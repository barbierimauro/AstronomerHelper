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
maxdecnorth=90+lat

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



vec_t1=[]
vec_t2=[]
s_ym=[]
query_dimm=[]
i=-1
for year in range(1994,2022):
	for month in range(1,13):
		i=i+1
		if(month<=9):
			t_start=str(year)+'-0'+str(month)  +'-01T00:00:00Z'
			s_ym.append(str(year)+'0'+str(month))
		else:
			t_start=str(year)+'-' +str(month)  +'-01T00:00:00Z'
			s_ym.append(str(year)+str(month))
		if(month<9):
			t_end=str(year)+'-0'+str(month+1)  +'-01T00:00:00Z'
		else:
			t_end=str(year)+'-' +str(month+1)  +'-01T00:00:00Z'
		if(month==12):
			t_end  =str(year+1)+'-01-01T00:00:00Z'
		#print(t_start,t_end)
		vec_t1.append(t_start)
		vec_t2.append(t_end)
		query_dimm_single="""SELECT start_date, ra, dec, airmass, integration, fwhm FROM asm.ambient_lasilla WHERE start_date > '"""+vec_t1[i]+"""' AND start_date< '"""+vec_t2[i]+"""'"""
		#print(query_dimm_single)
		query_dimm.append(query_dimm_single)

nmonths=len(vec_t1)
#print(nmonths)






MAXRECORDS=15000000 # maximum value of maxrec
SMALLREC=500
FULLQUERY=True


##nmonths=1
for nm_ in range(nmonths):
	print('-----------------------------------------------------------')
	print(vec_t1[nm_])
	print()
	t1 = timelib.time()
	if(FULLQUERY):
		rawframes = tapobs.search(query=query_dimm[nm_],maxrec=MAXRECORDS)
	else:
		rawframes = tapobs.search(query=query_dimm[nm_],maxrec=SMALLREC)
	t2 = timelib.time()
	print('query executed in ', t2-t1, 'sec')
	tab_dimm=rawframes.to_table()
	ndimm=len(tab_dimm)
	if(ndimm==0):
		continue
	else:
		print()
		print("Number of points:",ndimm)
		print()
		print('remove the last character Z from each string in the meteo dates ...')
		t1=timelib.time()
		t_ = [tmp_t[:-1] for tmp_t in tab_dimm['start_date']]
		tempo_dimm=Time(t_,format='isot',location=loc,scale='utc')
		t2=timelib.time()
		print('exec time ',t2-t1,'s')
		print()
		print('define variables from dimm ...')
		t1=timelib.time()
		airmass=np.array(tab_dimm['airmass'])
		ra=np.array(tab_dimm['ra'])
		dec=np.array(tab_dimm['dec'])
		fwhm=np.array(tab_dimm['fwhm'])
		t2=timelib.time()
		del tab_meteo
		print('exec time ',t2-t1,'s')
		print()
		t1=timelib.time()
		print('JD')
		jd=tempo_meteo.jd
		t2=timelib.time()
		print('exec time ',t2-t1,'s')
		print()

		object="COROT-11"
		#res=Simbad.query_object(object)
		#print(res )

		#creo le coordinate di un oggetto arbitrario
		object_c = coord.SkyCoord(0, 0, unit=(u.deg, u.deg), frame='icrs')

		glon=np.zeros(ndimm)
		glat=np.zeros(ndimm)
		elon=np.zeros(ndimm)
		elat=np.zeros(ndimm)
		az=np.zeros(ndimm)
		alt=np.zeros(ndimm)

		for i in range(ndimm):
			if(ra[i]>=0 and ra[i]<360 and dec[i]>=-90 and dec[i]<maxdecnorth):
				object_c[i] = coord.SkyCoord(ra[i], dec[i], unit=(u.deg, u.deg), frame='icrs')
				glon[i]=object_c[i].galactic.l.deg
				glat[i]=object_c[i].galactic.b.deg
				elon[i]=object_c[i].transform_to('barycentricmeanecliptic').lon.deg
				elat[i]=object_c[i].transform_to('barycentricmeanecliptic').lat.deg
				alt_az_frame = AltAz(location=loc, obstime=tempo_dimm[i])
				object_alt_az = object_c[i].transform_to(alt_az_frame)
				az[i] =object_alt_az.az.deg
				alt[i]=object_alt_az.alt.deg
				pointing[i] = Simbad.query_region(object_c[i],radius=0.17 * u.deg)
			else:
				alt[i]=90
				zenith_altaz = coord.SkyCoord(AltAz(az=0, alt=alt[i], obstime=tempo_dimm[i], location=loc))
				zenith_radec=zenith_altaz.transform_to(coord.ICRS)
				zenith_eclip=zenith_radec.transform_to('barycentricmeanecliptic')
				zenith_galac=zenith_radec.transform_to('galactic')
				ra[i]=zenith_radec.ra.deg
				dec[i]=zenith_radec.dec.deg
				glon[i]=zenith_galac.l.deg
				glat[i]=zenith_galac.b.deg
				elon[i]=zenith_eclip.lon.deg
				elat[i]=zenith_eclip.lat.deg
				object_c[i] = coord.SkyCoord(ra[i], dec[i], unit=(u.deg, u.deg), frame='icrs')
				pointing[i] = Simbad.query_region(object_c[i],radius=0.17 * u.deg)





		t1=timelib.time()
		print('writing')


		scrivi=1
		if (scrivi==1):
			t1 = timelib.time()
			ofile='meteo_ls.'+s_ym[nm_]+'.csv.gz'
			print(ofile)
			f_out=gzip.open(ofile,'wt',newline='')
			writer = csv.writer(f_out,delimiter=',')
			writer.writerow(['time','dt','JD',
							'year_decimal','year_frac','day_frac','day_frac_solar','year','month','dom<','dow','hour_ut','hour_solar','hour_night',
							'ra','dec','glon','glat','elon','elat','airmass','fwhm'])
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
				mag_sky[i],mag_sky_wmoon[i],flux_skybg_solaract[i],flux_skybg_zodiacal[i],flux_skybg_galactic[i],flux_skybg_moon[i],
				zenith_galac.l.value[i],zenith_galac.b.value[i],
				zenith_eclip.lon.value[i],zenith_eclip.lat.value[i]
				])

			f_out.close()

		t2=timelib.time()
		print('exec time ',t2-t1,'s')
		print()

#idx=np.where(daylight<=1)
#print(np.sum(idx)/nmeteo)
