#!/usr/bin/python3

##https://www.machinelearningplus.com/python/parallel-processing-python/

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
import multiprocessing as mp


# misc
as4rad=180.0*3600.0/np.pi
d2r=np.pi/180.0
r2d=180.0/np.pi




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




#La Silla
lon=-70.7375
lat=-29.2575
height=2400.0
loc=coord.EarthLocation(lat = lat*u.deg, lon = lon*u.deg, height=height*u.m)
delta_lon=(loc.lon.value/15)/24

cosmic_bg_mag_V=27.78 #https://www.ing.iac.es//Astronomy/observing/conditions/skybr/skybr.html

nfilters=9
filter_label=[    'U',   'B',   'V',   'R',   'I',   'Z',   'J',   'H',   'K']
filter_wav  =[   3600,  4300,  5500,  6500,  8200,  9500, 12500, 16700, 21600]
filter_band =[    500,   720,   860,  1330,  1400,   700,  1600,  2900,  3200]
sky_flux_jy =[   1810,  4260,  3640,  3080,  2550,  2200,  1570,  1020,   640]
sky_flux_phs=[    783,  1463,  1050,   873,   376,   329,   194,  94.4,  44.1]
atmo_extin  =[   0.55,  0.25,  0.15,  0.09,  0.06,  0.05,  0.10,  0.11,  0.07]
skymagbg=[       22.0,  22.7,  21.9,  21.0,  20.0,  18.8,  16.1,  14.7,  12.5]
color_label =[  'U-B', 'B-V', 'V-R', 'R-I', 'I-Z', 'Z-J', 'J-H', 'H-K','NONE']
skymagcolor =[   -0.7,   0.8,   0.9,   1.0,   1.2,   2.7,   1.4,   2.2,np.nan]
fudge=[           2.0,   1.6,   2.4,   3.0,   5.9,   0.0,   0.0,   0.0,   0.0]




MAXRECORDS=15000000 # maximum value of maxrec
FULLQUERY=False
SMALLREC=5000

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
del tab_meteo
print('exec time ',t2-t1,'s')
print()




t1=timelib.time()
print('JD,LST')
jd=tempo_meteo.jd
##pool = mp.Pool(mp.cpu_count())
##jd=[pool.apply(tempo_meteo.jd), []]
##pool.close()
##sys.exit()
#mjd=tempo_meteo.mjd
lst=tempo_meteo.sidereal_time('apparent')
t2=timelib.time()
print('exec time ',t2-t1,'s')
print()

t1=timelib.time()
print('AltAz,RaDec,Gal,Ecl')
# converte az e alt in RA,DEC
# non si puo' fare come si fa con sole e luna perche' AltAz e' buggato e non accetta vettori di coordinate ma solo una unica coordinata
# si risolve creando dei vettori con le coordinate
# in questo caso creo az e alt dello zenith e li converto in RA DEC
az = np.zeros(nmeteo)
alt = np.ones(nmeteo) * 90.0
# definisco il frame delle coordinate
zenith_altaz = coord.SkyCoord(AltAz(az=az * u.deg, alt=alt * u.deg, obstime=tempo_meteo, location=loc))
# converto le coordinate in RA,DEC
# attenzione che coord.ICRS e' il package coordinate di astropy e ICRS e' uno dei frame che gestisce
# lista dei frames
#https://docs.astropy.org/en/stable/coordinates/index.html#built-in-frame-classes

zenith_radec=zenith_altaz.transform_to(coord.ICRS)
zenith_eclip=zenith_radec.transform_to('barycentricmeanecliptic')
zenith_galac=zenith_radec.transform_to('galactic')
del zenith_altaz
del alt
del az
del zenith_radec
t2=timelib.time()
print('exec time ',t2-t1,'s')
print()

#t1=timelib.time()
#print('BJD')
#ltt_bary = tempo_meteo.light_travel_time(zenith_radec).tdb.jd
#bjd_tdb=tempo_meteo.tdb.jd + ltt_bary
#ltt_bary = tempo_meteo.light_travel_time(zenith_radec).tcb.jd
#bjd_tcb=tempo_meteo.tcb.jd + ltt_bary
#ltt_bary=0
#t2=timelib.time()
#print('exec time ',t2-t1,'s')
#print()

t1=timelib.time()
print('Sun')
alt_az_frame = AltAz(location=loc, obstime=tempo_meteo)
sun=get_sun(tempo_meteo)
sun_coord=coord.SkyCoord(sun.ra,sun.dec, unit=(u.deg,u.deg))
sun_alt_az = sun_coord.transform_to(alt_az_frame)
#sep_sun=zenith_radec.separation(sun_coord)
dist_sun=sun.distance
t2=timelib.time()
print('exec time ',t2-t1,'s')
print()

t1=timelib.time()
print('Moon')
moon=get_moon(tempo_meteo,loc)
print('moon coord')
moon_coord=coord.SkyCoord(moon.ra,moon.dec, unit=(u.deg,u.deg))
print('moon elongation')
moon_elongation = sun.separation(moon_coord)
print('dist moon')
dist_moon=moon.distance
print('moon d1')
d1=dist_sun.to(u.km)*np.sin(moon_elongation.rad)
print('moon d2')
d2=dist_moon - d1*np.cos(moon_elongation.rad)
print('moon phase rad')
moon_phase_rad=np.arctan2(d1,d2)
del d1
del d2
print('moon phase angle')
moon_phase_angle=moon_phase_rad.value*r2d # 180=NM, 90=Q, 0=FM
print('moon altaz')
moon_alt_az = moon_coord.transform_to(alt_az_frame)
print('moon sep')
sep_moon=90-moon_alt_az.alt.deg

# clear useless vectors
del sun
del sun_coord
del moon
del moon_coord
del moon_phase_rad
del moon_elongation

t2=timelib.time()
print('exec time ',t2-t1,'s')
print()


t1=timelib.time()
print('daylight,moonlight')
daylight=np.zeros(nmeteo)
moonlight=np.zeros(nmeteo)
# astronomical twilight
idx=np.where((sun_alt_az.alt.deg>=-18) & (sun_alt_az.alt.deg<-12))
daylight[idx]=1
# nautical twilight
idx=np.where((sun_alt_az.alt.deg>=-12) & (sun_alt_az.alt.deg<-6))
daylight[idx]=2
# civil twilight
idx=np.where((sun_alt_az.alt.deg>=-6) & (sun_alt_az.alt.deg<-1))
daylight[idx]=3
# sun above horizon
idx=np.where(sun_alt_az.alt.deg>=-1)
daylight[idx]=4

# moon twilight
idx=((moon_alt_az.alt.deg>=-12) & (moon_alt_az.alt.deg<-1))
moonlight[idx]=1
# moon above horiz
idx=(moon_alt_az.alt.deg>=-1)
moonlight[idx]=2


t2=timelib.time()
print('exec time ',t2-t1,'s')
print()





t1=timelib.time()
print('sky background')

solar_cycle_phase=np.zeros(nmeteo)
# calculate solar phase
#solar minima: 1986.75, 1996.67, 2009.0, 2020.0
#period lenght             9.92,  12.33,   11.0
#
# cycle 21
idx=np.where((tempo_meteo.jyear>1986.75) & (tempo_meteo.jyear<= 1996.67))
cycle_len=9.92
solar_cycle_phase[idx] = (tempo_meteo.jyear[idx]-1986.75)/cycle_len
# cycle 22
idx=np.where((tempo_meteo.jyear>1996.67) & (tempo_meteo.jyear<= 2009))
cycle_len=12.33
solar_cycle_phase[idx] = (tempo_meteo.jyear[idx]-1996.67)/cycle_len
# cycle 23
idx=np.where((tempo_meteo.jyear>2009) & (tempo_meteo.jyear<= 2020))
cycle_len=11
solar_cycle_phase[idx] = (tempo_meteo.jyear[idx]-2009)/cycle_len
# cycle 24
idx=np.where(tempo_meteo.jyear>2020)
cycle_len=11
solar_cycle_phase[idx] = (tempo_meteo.jyear[idx]-2020)/cycle_len


flux_skybg_solaract=np.zeros(nmeteo)
flux_skybg_zodiacal=np.zeros(nmeteo)
flux_skybg_galactic=np.zeros(nmeteo)
flux_skybg_moon    =np.zeros(nmeteo)
flux_skybg_total   =np.zeros(nmeteo)
mag_sky            =np.full(nmeteo,np.nan)
mag_sky_wmoon      =np.full(nmeteo,np.nan)

# during day set the sky magnitude
idx=np.where(daylight>1)
mag_sky[idx]      =-8
mag_sky_wmoon[idx]=-8
flux_skybg_solaract[idx]=np.nan
flux_skybg_zodiacal[idx]=np.nan
flux_skybg_galactic[idx]=np.nan
flux_skybg_moon[idx]    =np.nan

# when is night calculate
# solar activity background and
# milky way light
idx=np.where(daylight<=1)
flux_skybg_solaract[idx]=(145+130*(solar_cycle_phase[idx]-0.8)/1.2)
#skybg_solaract=58.3 # at solar minimum
flux_skybg_galactic=100*np.exp(-np.abs(zenith_galac.b.value)/10)

#when is night and the ecliptic latitude is below 60 deg calculate the
#zodiacal light
idx=np.where((daylight<=1)&(zenith_eclip.lat.value>=60))
flux_skybg_zodiacal[idx]=60
idx=np.where((daylight<=1)&(zenith_eclip.lat.value<60))
flux_skybg_zodiacal[idx]=140-90*np.sin(zenith_eclip.lat.rad[idx])

#calculate dark-sky brightness
idx=np.where(daylight<=1)
mag_sky[idx]=cosmic_bg_mag_V-2.5*np.log10(flux_skybg_solaract[idx] + flux_skybg_zodiacal[idx] + flux_skybg_galactic[idx])

# during the night when the moon is below the horizon, the sky is already calculated
idx=np.where((daylight<=1)&(moonlight==0))
mag_sky_wmoon[idx]=mag_sky[idx]

# moon above horizon during the night
idx=np.where((daylight<=1)&(moonlight>0))
s   =np.zeros(nmeteo)
fr  =np.zeros(nmeteo)
xzm =np.zeros(nmeteo)
bs10=np.zeros(nmeteo)

s[idx]  = 10**(-0.4*(3.84+0.026*moon_phase_angle[idx]+4e-9*moon_phase_angle[idx]**4))
fr[idx] = 10**5.36*(1.06+(np.cos(sep_moon[idx]*d2r))**2) +10**(6.15-sep_moon[idx]/40)
#xzm[idx]= 1/np.sqrt(1-0.96*(np.sin((90-moon_alt_az.alt.value[idx])*d2r))**2)
xzm[idx]= 1/np.sqrt(1-0.96*(np.cos(moon_alt_az.alt.rad[idx]))**2)
ak      = atmo_extin[2]
k1      = 3.8*(1-10**(-0.4*ak))*fudge[2]
flux_skybg_moon[idx] = k1*s[idx]*fr[idx]*10**(-0.4*ak*xzm[idx])
mag_sky_wmoon[idx] = cosmic_bg_mag_V-2.5*np.log10(flux_skybg_moon[idx] + flux_skybg_solaract[idx] + flux_skybg_zodiacal[idx] + flux_skybg_galactic[idx])

del s
del fr
del xzm
del sep_moon

t2=timelib.time()
print('exec time ',t2-t1,'s')
print()





t1=timelib.time()
print('writing')


scrivi=1
if (scrivi==1):
	t1 = timelib.time()
	f_out=gzip.open('meteo_ls.csv.gz','wt',newline='')
	writer = csv.writer(f_out,delimiter=',')
	writer.writerow(['time','dt','JD','LST',
					'year_decimal','year_frac','day_frac','day_frac_solar','year','month','dom','dow','hour_ut','hour_solar','hour_night',
					'pressure2m','relhum2m','temp2m','tempdew2m','winddir10m','windspeed10m','temp30m','winddir30m','windspeed30m',
					'sun_alt','sun_az','sun_dist','moon_alt','moon_az','moon_phase','moon_dist',
					'daylight','moonlight','solar_cycle_phase',
					'mag_sky','mag_sky_wmoon','skybg_solaract','skybg_zodiacal','skybg_galactic','skybg_moon',
					'zenith_l','zenith_b','zenith_elon','zenith_elat'])
	for i in range(nmeteo):
		d_t           =tempo_meteo.jd[i]-tempo_meteo.jd[i-1]
		if(i==0):
			d_t=0
		year_decimal  =tempo_meteo.jyear[i]
		year_frac     =np.mod(tempo_meteo.jyear[i],1)
		day_frac      =np.mod(tempo_meteo.jd[i]+0.5,1)
		day_frac_solar=np.mod(day_frac+delta_lon,1)
		hour_ut       =int(tempo_meteo.isot[i][11:13])
		hour_solar    =int(day_frac_solar*24)
		hour_night    =(np.mod(day_frac_solar+0.5,1)-0.5)*24
		year          =int(tempo_meteo.isot[i][0:4])
		month         =int(tempo_meteo.isot[i][5:7])
		dom           =int(tempo_meteo.isot[i][8:10])
		dow           =datetime.datetime(year,month,dom).weekday()

		writer.writerow([
		tempo_meteo[i],d_t,jd[i],lst.value[i],
		year_decimal,year_frac,day_frac,day_frac_solar,year,month,dom,dow,hour_ut,hour_solar,hour_night,
		pressure2m[i],rhum2m[i],temp2m[i],tempdew2m[i],winddir10m[i],windspeed10m[i],
		temp30m[i],winddir30m[i],windspeed30m[i],
		sun_alt_az.alt.value[i],sun_alt_az.az.value[i],dist_sun.value[i],
		moon_alt_az.alt.value[i],moon_alt_az.az.value[i],(180-moon_phase_angle[i])/180,dist_moon.value[i],
		daylight[i],moonlight[i],solar_cycle_phase[i],
		mag_sky[i],mag_sky_wmoon[i],flux_skybg_solaract[i],flux_skybg_zodiacal[idx],flux_skybg_galactic[i],flux_skybg_moon[i],
		zenith_galac.l.value[i],zenith_galac.b.value[i],
		zenith_eclip.lon.value[i],zenith_eclip.lat.value[i]
		])

	f_out.close()

t2=timelib.time()
print('exec time ',t2-t1,'s')
print()

idx=np.where(daylight<=1)
print(np.sum(idx)/nmeteo)
