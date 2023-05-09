#!/usr/bin/python3


import os
import sys
import math
import pyvo
from pyvo.dal import tap
import requests
import cgi
import re
import json
import getpass
import numpy as np
import scipy.stats as st
from statsmodels.sandbox.stats import runs
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
import numpy.ma as ma
import astropy.coordinates as coord
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation, AltAz
from astropy.coordinates import get_body_barycentric, get_body, get_moon, get_sun
from scipy import interpolate
import csv
import gzip
import time as timelib
import datetime

# misc
as4rad=180.0*3600.0/np.pi
d2r=np.pi/180.0
r2d=180.0/np.pi


# ESO VO definitions

mode_requested = "raw2master"
thisscriptname = "jupyter_script"
headers={'User-Agent': '%s (ESO script drc %s)'%(requests.utils.default_headers()['User-Agent'], thisscriptname)}
def getToken(username, password):
    """Token based authentication to ESO: provide username and password to receive back a JSON Web Token."""
    if username==None or password==None:
        return None
    token_url = "https://www.eso.org/sso/oidc/token"
    token = None
    try:
        response = requests.get(token_url,
                            params={"response_type": "id_token token", "grant_type": "password",
                                    "client_id": "clientid",
                                    "username": username, "password": password})
        token_response = json.loads(response.content)
        token = token_response['id_token']+'=='
    except NameError as e:
        print(e)
    except:
        print("*** AUTHENTICATION ERROR: Invalid credentials provided for username %s" %(username))

    return token

def download_asset(file_url, filename=None, token=None):

    headers = None
    if token!=None:
        headers = {"Authorization": "Bearer " + token}
        response = requests.get(file_url, headers=headers)
    else:
        # Trying to download anonymously
        response = requests.get(file_url, stream=True, headers=headers)

    if filename == None:
        contentdisposition = response.headers.get('Content-Disposition')
        if contentdisposition != None:
            value, params = cgi.parse_header(contentdisposition)
            filename = params["filename"]

        if filename == None:
            # last chance: get anything after the last '/'
            filename = url[url.rindex('/')+1:]

    if response.status_code == 200:
        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=50000):
                f.write(chunk)

    return (response.status_code, filename)


# ESO credentials and db access

username='mbarbierops'
password='ACenbb2022#'
token = getToken(username, password)
if token == None:
    print("Could not authenticate. Continuing as anonymous")
else:
    print("Authentication successful")
print()


ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
tapobs = tap.TAPService(ESO_TAP_OBS)

ESO_TAP_CAT = "http://archive.eso.org/tap_cat"
tapcat = tap.TAPService(ESO_TAP_CAT)



#astronomical site definitions

#La Silla
lon=-70.7375
lat=-29.2575
height=2400.0
loc=coord.EarthLocation(lat = lat*u.deg, lon = lon*u.deg, height=height*u.m)
delta_lon=(loc.lon.value/15)/24

darksky=27.78 #https://www.ing.iac.es//Astronomy/observing/conditions/skybr/skybr.html
light_pollution_mag=0.5

nfilters=9
skymag_filter_label=[    'U',   'B',   'V',   'R',   'I',   'Z',   'J',   'H',   'K']
skymag_filter_wav=[     3600,  4300,  5500,  6500,  8200,  9500, 12500, 16700, 21600]
skymag_filter_band=[     500,   720,   860,  1330,  1400,   700,  1600,  2900,  3200]
skymag_filter_flux_jy=[ 1810,  4260,  3640,  3080,  2550,  2200,  1570,  1020,   640]
skymag_filter_flux_phs=[ 783,  1463,  1050,   873,   376,   329,   194,  94.4,  44.1]
skymag_filter_extin=[   0.55,  0.25,  0.15,  0.09,  0.06,  0.05,  0.10,  0.11,  0.07]
skymag=[                22.0,  22.7,  21.9,  21.0,  20.0,  18.8,  16.1,  14.7,  12.5]
fudge=[                  2.0,   1.6,   2.4,   3.0,   5.9,   0.0,   0.0,   0.0,   0.0]


# Query

print('start query')
MAXRECORDS=15000000 # maximum value of maxrec
FULLQUERY=True
SMALLREC=1000

query_meteo="""
SELECT start_date,press_2m,rhum_2m,temp_2m,tempdew_2m,temp_30m,wind_dir_30m,wind_speed_30m,wind_dir_10m,wind_speed_10m
from asm.meteo_lasilla
where start_date > '1994-01-01T00:00:00Z' and start_date< '2022-01-01T00:00:00Z'
"""
t1 = timelib.time()
if(FULLQUERY):
    rawframes = tapobs.search(query=query_meteo,maxrec=MAXRECORDS)
else:
    rawframes = tapobs.search(query=query_meteo,maxrec=SMALLREC)
t2 = timelib.time()
print('query executed in ', t2-t1, 'sec')
tab_meteo=rawframes.to_table()
nmeteo=len(tab_meteo)
print()
print('----------------------')
print("Metereological dataset")
print('----------------------')
print("Number of points:",nmeteo)
print()
print("First data row")
print(tab_meteo[0])
print()



print('remove the last character Z from each string in the meteo dates ...')
t1=timelib.time()
t_ = [tmp_t[:-1] for tmp_t in tab_meteo['start_date']]
tempo_meteo=Time(t_,format='isot',location=loc,scale='utc')
t2=timelib.time()
print('exec time ',t2-t1,'s')
print()
print('define variables from meteo ...')
t1=timelib.time()
pressure2m  =np.array(tab_meteo['press_2m'])
rhum2m      =np.array(tab_meteo['rhum_2m'])
temp2m      =np.array(tab_meteo['temp_2m'])
tempdew2m   =np.array(tab_meteo['tempdew_2m'])
winddir10m  =np.array(tab_meteo['wind_dir_10m'])
windspeed10m=np.array(tab_meteo['wind_speed_10m'])
temp30m     =np.array(tab_meteo['temp_30m'])
winddir30m  =np.array(tab_meteo['wind_dir_30m'])
windspeed30m=np.array(tab_meteo['wind_speed_30m'])
t2=timelib.time()
print('exec time ',t2-t1,'s')
print()

del(tab_meteo)
del(t_)

f_out=gzip.open('meteo_ls.csv.gz','wt',newline='')
writer = csv.writer(f_out,delimiter=',')
writer.writerow(['time','dt','JD','MJD','BJD_TDB','BJD_TCB','LST',
				 'year_decimal','year_frac','day_frac','day_frac_solar','year','month','dom','dow','hour_ut','hour_solar','hour_night',
				 'pressure2m','relhum2m','temp2m','tempdew2m','winddir10m','windspeed10m','temp30m','winddir30m','windspeed30m',
				 'sun_alt','sun_az','sun_dist','moon_alt','moon_az','moon_phase','moon_dist','daylight','moonlight','solar_cycle_phase',
				 'sky_mag_U','sky_mag_B','sky_mag_V','sky_mag_R','sky_mag_I','sky_mag_Z','sky_mag_J','sky_mag_H','sky_mag_K',
				 'sky_mag_U_Moon','sky_mag_B_Moon','sky_mag_V_Moon','sky_mag_R_Moon','sky_mag_I_Moon','sky_mag_Z_Moon','sky_mag_J_Moon','sky_mag_H_Moon','sky_mag_K_Moon',
				 'skybg_airm','skybg_zod','skybg_mw','skybg_moon','zenith_ra','zenith_dec','zenith_l','zenith_b','zenith_elon','zenith_elat'])

t1=timelib.time()
print('writing')


for i in range(nmeteo):
	if(i % 10 ==0):
		print(i, nmeteo, i/nmeteo*100,'%')
	jd =tempo_meteo.jd[i]
	mjd=tempo_meteo.mjd[i]
	lst=tempo_meteo.sidereal_time('apparent')[i]
	zenith_altaz = coord.SkyCoord(AltAz(az=0.0*u.deg, alt=90.0*u.deg, obstime=tempo_meteo[i], location=loc))
	zenith_radec = zenith_altaz.transform_to(coord.ICRS)
	zenith_eclip = zenith_radec.transform_to('barycentricmeanecliptic')
	zenith_galac = zenith_radec.transform_to('galactic')
	ltt_bary     = tempo_meteo.light_travel_time(zenith_radec).tdb.jd[i]
	bjd_tdb      = tempo_meteo.tdb.jd[i] + ltt_bary
	ltt_bary     = tempo_meteo.light_travel_time(zenith_radec).tcb.jd[i]
	bjd_tcb      = tempo_meteo.tcb.jd[i] + ltt_bary
	alt_az_frame = AltAz(location=loc, obstime=tempo_meteo[i])
	sun=get_sun(tempo_meteo[i])
	sun_coord=coord.SkyCoord(sun.ra,sun.dec, unit=(u.deg,u.deg))
	sun_alt_az = sun_coord.transform_to(alt_az_frame)
	dist_sun=sun.distance
	moon=get_moon(tempo_meteo[i],loc)
	moon_coord=coord.SkyCoord(moon.ra,moon.dec, unit=(u.deg,u.deg))
	moon_elongation = sun.separation(moon_coord)
	dist_moon=moon.distance
	moon_phase_rad=np.arctan2(dist_sun.to(u.km)*np.sin(moon_elongation.rad),dist_moon - dist_sun.to(u.km)*np.cos(moon_elongation.rad))
	moon_phase_angle=moon_phase_rad.value*r2d # 180=NM, 90=Q, 0=FM
	moon_alt_az = moon_coord.transform_to(alt_az_frame)
	sep_moon=zenith_radec.separation(moon_coord)
	alts=sun_alt_az.alt.deg
	if(alts<-18):
		daylight=0
	elif(alts>=-18 and alts<-12):
		daylight=1
	elif(alts>=-12 and alts<-6):
		daylight=2
	elif(alts>=-6):
		daylight=3
	altm=moon_alt_az.alt.deg
	if(altm<0):
		moonlight=0
	elif(altm>=-6 and altm<0):
		moonlight=1
	elif(altm>=0):
		moonlight=2
#solar minima: 1986.75, 1996.67, 2009.0, 2020.0
#period lenght             9.92,  12.33,   11.0
	if(tempo_meteo.jyear[i]>1986.75 and tempo_meteo.jyear[i]<= 1996.67):
		cycle_len=9.92
		solar_cycle_phase = (tempo_meteo.jyear[i]-1986.75)/cycle_len
	elif(tempo_meteo.jyear[i]>1996.67 and tempo_meteo.jyear[i]<= 2009.0):
		cycle_len=12.33
		solar_cycle_phase = (tempo_meteo.jyear[i]-1996.67)/cycle_len
	elif(tempo_meteo.jyear[i]>2009 and tempo_meteo.jyear[i]<= 2020):
		cycle_len=11
		solar_cycle_phase = (tempo_meteo.jyear[i]-2009)/cycle_len
	elif(tempo_meteo.jyear[i]>2020):
		cycle_len=11
		solar_cycle_phase = (tempo_meteo.jyear[i]-2020)/cycle_len

	sky_mag      =np.full(nmeteo,np.nan)
	sky_mag_wmoon=np.full(nmeteo,np.nan)

	ak=skymag_filter_extin[2]
	xz=1

	if(daylight>1):
		sky_mag[:]      =-8
		sky_mag_wmoon[:]=-8
	else:
		# solar activity background
		skybg_airm=(145+130*(solar_cycle_phase-0.8)/1.2)
		#skybg_airm=58.3 # at solar minimum
		# zodiacal light
		if(zenith_eclip.lat.value<60):
			skybg_zod=140-90*np.sin(zenith_eclip.lat.value*d2r)
		else:
			skybg_zod=60
		#skybg_zod=60
		# milky way light
		skybg_mw=100*np.exp(-np.abs(zenith_galac.b.value)/10)
		skybg_tot=skybg_airm + skybg_zod + skybg_mw
		#   Assume sky colour implied by values in SKYMAG, and calculate
		#   dark-sky brightness in the other bands

		skylocal=np.full(nfilters,np.nan)
		qskylocal=np.zeros(nfilters)
		for j in range(nfilters):
			skylocal[j]=darksky-2.5*np.log10(skybg_tot)+skymag[j]-skymag[2]-light_pollution_mag
			qskylocal[j]=10.**((darksky-skylocal[j])/2.5)
			sky_mag[j]=skylocal[j]

		skylocalwmoom=np.full(nfilters,np.nan)
		qmoon=np.zeros(nfilters)
		qall=np.zeros(nfilters)
		###zm=90-moon_alt_az.alt.value
		if(moonlight==0):
			for j in range(nfilters):
				skylocalwmoom[j]=skylocal[j]
		else:
			#xz =airmass
			s=10**(-0.4*(3.84+0.026*moon_phase_angle+4e-9*moon_phase_angle**4))
			fr=10**5.36*(1.06+(np.cos(sep_moon.value*d2r))**2)+10**(6.15-sep_moon.value/40)
			xzm=1/np.sqrt(1-0.96*(np.sin((90-moon_alt_az.alt.value)*d2r))**2)
			bnl=s*fr*10**(-0.4*ak*xzm)*(1-10**(-0.4*ak*xz))
			bs10=bnl*3.8
			for j in range(nfilters):
				qmoon[j] =fudge[j]*bs10
				qall[j]  =qskylocal[j]+qmoon[j]
				skylocalwmoom[j]=darksky-2.5*np.log10(qall[j])+skymag[j]-skymag[2]-light_pollution_mag
				sky_mag_wmoon[j]=skylocalwmoom[j]

		skybg_moon=qmoon[2]


		if(i==0):
			d_t=0
		else:
			d_t=(tempo_meteo.jd[i]-tempo_meteo.jd[i-1])*86400
		year_decimal  =tempo_meteo.jyear[i]
		year_frac     =np.mod(tempo_meteo.jyear[i],1)
		day_frac      =np.mod(tempo_meteo.jd[i]+0.5,1)
		day_frac_solar=np.mod(day_frac+delta_lon,1)
		hour_ut       =int(tempo_meteo.isot[i][11:13])
		hour_solar    =int(day_frac_solar*24)
		year          =int(tempo_meteo.isot[i][0:4])
		month         =int(tempo_meteo.isot[i][5:7])
		dom           =int(tempo_meteo.isot[i][8:10])
		dow           =datetime.datetime(year,month,dom).weekday()
		hour_night    =(np.mod(day_frac_solar+0.5,1)-0.5)*24

		writer.writerow([
		tempo_meteo[i],d_t,jd,mjd,bjd_tdb,bjd_tcb,lst.value,
		year_decimal,year_frac,day_frac,day_frac_solar,year,month,dom,dow,hour_ut,hour_solar,hour_night,
		pressure2m[i],rhum2m[i],temp2m[i],tempdew2m[i],winddir10m[i],windspeed10m[i],
		temp30m[i],winddir30m[i],windspeed30m[i],
		sun_alt_az.alt.value,sun_alt_az.az.value,dist_sun.value,
		moon_alt_az.alt.value,moon_alt_az.az.value,(180-moon_phase_angle)/180,dist_moon.value,
		daylight,moonlight,solar_cycle_phase,
		sky_mag[0],sky_mag[1],sky_mag[2],sky_mag[3],sky_mag[4],sky_mag[5],sky_mag[6],sky_mag[7],sky_mag[8],
		sky_mag_wmoon[0],sky_mag_wmoon[1],sky_mag_wmoon[2],sky_mag_wmoon[3],sky_mag_wmoon[4],sky_mag_wmoon[5],sky_mag_wmoon[6],sky_mag_wmoon[7],sky_mag_wmoon[8],
		skybg_airm,skybg_zod,skybg_mw,skybg_moon,
		zenith_radec.ra.value,zenith_radec.dec.value,
		zenith_galac.l.value,zenith_galac.b.value,
		zenith_eclip.lon.value,zenith_eclip.lat.value
		])















f_out.close()
t2=timelib.time()
print('exec time ',t2-t1,'s')
print()
