#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


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


# In[3]:


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


# In[4]:


#La Silla
lon=-70.7375
lat=-29.2575
height=2400.0
loc=coord.EarthLocation(lat = lat*u.deg, lon = lon*u.deg, height=height*u.m)

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


# In[5]:


MAXRECORDS=15000000 # maximum value of maxrec
FULLQUERY=True
SMALLREC=2000

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



# In[6]:


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


# In[7]:


t1=timelib.time()
print('JD,MJD,LST')
jd=tempo_meteo.jd
mjd=tempo_meteo.mjd
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
t2=timelib.time()
print('exec time ',t2-t1,'s')
print()

t1=timelib.time()
print('BJD')
ltt_bary = tempo_meteo.light_travel_time(zenith_radec).tdb.jd
bjd_tdb=tempo_meteo.tdb.jd + ltt_bary
ltt_bary = tempo_meteo.light_travel_time(zenith_radec).tcb.jd
bjd_tcb=tempo_meteo.tcb.jd + ltt_bary
ltt_bary=0
t2=timelib.time()
print('exec time ',t2-t1,'s')
print()

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
moon_coord=coord.SkyCoord(moon.ra,moon.dec, unit=(u.deg,u.deg))
moon_elongation = sun.separation(moon_coord)
dist_moon=moon.distance
moon_phase_rad=np.arctan2(dist_sun.to(u.km)*np.sin(moon_elongation.rad),dist_moon - dist_sun.to(u.km)*np.cos(moon_elongation.rad))
moon_phase_angle=moon_phase_rad.value*r2d # 180=NM, 90=Q, 0=FM
moon_alt_az = moon_coord.transform_to(alt_az_frame)
sep_moon=zenith_radec.separation(moon_coord)
##print(sun.distance[0].to(u.km))
##print(moon.distance[0])

# clear useless vectors
az=0
alt=0
sun=0
sun_coord=0
moon=0
moon_coord=0
moon_phase_rad=0
moon_elongation=0

t2=timelib.time()
print('exec time ',t2-t1,'s')
print()


# In[8]:



t1=timelib.time()
print('daylight,moonlight')
daylight=np.zeros(nmeteo)
moonlight=np.zeros(nmeteo)
idx=np.where((sun_alt_az.alt.deg>=-18) & (sun_alt_az.alt.deg<-12))
daylight[idx]=1
idx=np.where((sun_alt_az.alt.deg>=-12) & (sun_alt_az.alt.deg<-6))
daylight[idx]=2
idx=np.where((sun_alt_az.alt.deg>=-6) & (sun_alt_az.alt.deg<-1))
daylight[idx]=3
idx=np.where(sun_alt_az.alt.deg>=-1)
daylight[idx]=4

idx=((moon_alt_az.alt.deg>=-12) & (moon_alt_az.alt.deg<-1))
moonlight[idx]=1
idx=(moon_alt_az.alt.deg>=-1)
moonlight[idx]=2


t2=timelib.time()
print('exec time ',t2-t1,'s')
print()




# In[9]:


t1=timelib.time()
print('sky background')

skybg_airm=np.full(nmeteo,np.nan)
skybg_zod =np.full(nmeteo,np.nan)
skybg_mw  =np.full(nmeteo,np.nan)
skybg_moon=np.full(nmeteo,np.nan)
skybg_tot =np.full(nmeteo,np.nan)
sky_mag      =np.full((nmeteo,nfilters),np.nan)
sky_mag_wmoon=np.full((nmeteo,nfilters),np.nan)
solar_cycle_phase=np.zeros(nmeteo)

ak=skymag_filter_extin[2]
xz=1

for i in range(nmeteo):

#solar minima: 1986.75, 1996.67, 2009.0, 2020.0
#period lenght             9.92,  12.33,   11.0
    if(tempo_meteo.jyear[i]>1986.75 and tempo_meteo.jyear[i]<= 1996.67):
        cycle_len=9.92
        solar_cycle_phase[i] = (tempo_meteo.jyear[i]-1986.75)/cycle_len
    elif(tempo_meteo.jyear[i]>1996.67 and tempo_meteo.jyear[i]<= 2009.0):
        cycle_len=12.33
        solar_cycle_phase[i] = (tempo_meteo.jyear[i]-1996.67)/cycle_len
    elif(tempo_meteo.jyear[i]>2009 and tempo_meteo.jyear[i]<= 2020):
        cycle_len=11
        solar_cycle_phase[i] = (tempo_meteo.jyear[i]-2009)/cycle_len
    elif(tempo_meteo.jyear[i]>2020):
        cycle_len=11
        solar_cycle_phase[i] = (tempo_meteo.jyear[i]-2020)/cycle_len


    if(daylight[i]>1):
        sky_mag[i][:]      =-8
        sky_mag_wmoon[i][:]=-8

    else:
        # solar activity background
        skybg_airm[i]=(145+130*(solar_cycle_phase[i]-0.8)/1.2)
        #skybg_airm[i]=58.3 # at solar minimum

        # zodiacal light
        if(zenith_eclip.lat.value[i]<60):
            skybg_zod[i]=140-90*np.sin(zenith_eclip.lat.value[i]*d2r)
        else:
            skybg_zod[i]=60
        #skybg_zod[i]=60

        # milky way light
        skybg_mw[i]=100*np.exp(-np.abs(zenith_galac.b.value[i])/10)

        skybg_tot[i]=skybg_airm[i] + skybg_zod[i] + skybg_mw[i]


        #   Assume sky colour implied by values in SKYMAG, and calculate
        #   dark-sky brightness in the other bands

        skylocal=np.full(nfilters,np.nan)
        qskylocal=np.zeros(nfilters)
        for j in range(nfilters):
            skylocal[j]=darksky-2.5*np.log10(skybg_tot[i])+skymag[j]-skymag[2]-light_pollution_mag
            qskylocal[j]=10.**((darksky-skylocal[j])/2.5)
            sky_mag[i][j]=skylocal[j]


        #print(i,skybg_tot[i],skybg_airm[i],skybg_zod[i],skybg_mw[i],sky[2],solar_cycle_phase)
        #print(i,sky[2])



        skylocalwmoom=np.full(nfilters,np.nan)
        qmoon=np.zeros(nfilters)
        qall=np.zeros(nfilters)
        ###zm=90-moon_alt_az.alt.value
        if(moonlight[i]==0):
            for j in range(nfilters):
                skylocalwmoom[j]=skylocal[j]
        else:
            #xz =airmass
            s=10**(-0.4*(3.84+0.026*moon_phase_angle[i]+4e-9*moon_phase_angle[i]**4))
            fr=10**5.36*(1.06+(np.cos(sep_moon.value[i]*d2r))**2)+10**(6.15-sep_moon.value[i]/40)
            xzm=1/np.sqrt(1-0.96*(np.sin((90-moon_alt_az.alt.value[i])*d2r))**2)
            bnl=s*fr*10**(-0.4*ak*xzm)*(1-10**(-0.4*ak*xz))
            bs10=bnl*3.8
            for j in range(nfilters):
                qmoon[j] =fudge[j]*bs10
                qall[j]  =qskylocal[j]+qmoon[j]
                skylocalwmoom[j]=darksky-2.5*np.log10(qall[j])+skymag[j]-skymag[2]-light_pollution_mag
                sky_mag_wmoon[i][j]=skylocalwmoom[j]
                #print(i,moon_alt_az.alt.value[i],sky[2],skyall[2])

        skybg_moon[i]=qmoon[2]

skybg_tot=0
sep_moon=0

t2=timelib.time()
print('exec time ',t2-t1,'s')
print()



# In[10]:


t1=timelib.time()
print('writing')

delta_lon=(loc.lon.value/15)/24

scrivi=1
if (scrivi==1):
    t1 = timelib.time()
    f_out=gzip.open('meteo_ls.csv.gz','wt',newline='')
    writer = csv.writer(f_out,delimiter=',')
    writer.writerow(['time','dt','JD','MJD','BJD_TDB','BJD_TCB','LST',
                     'year_decimal','year_frac','day_frac','day_frac_solar','year','month','dom','dow','hour_ut','hour_solar','hour_night',
                     'pressure2m','relhum2m','temp2m','tempdew2m','winddir10m','windspeed10m','temp30m','winddir30m','windspeed30m',
                     'median_pressure2m','median_rhum2m','median_temp2m','median_tempdew2m','median_winddir10m','median_windspeed10m','median_temp30m','median_winddir30m','median_windspeed30m',
                     'd_pressure2m','d_rhum2m','d_temp2m','d_tempdew2m','d_winddir10m','d_windspeed10m','d_temp30m','d_winddir30m','d_windspeed30m',
                     'sun_alt','sun_az','sun_dist','moon_alt','moon_az','moon_phase','moon_dist','daylight','moonlight','solar_cycle_phase',
                     'sky_mag_U','sky_mag_B','sky_mag_V','sky_mag_R','sky_mag_I','sky_mag_Z','sky_mag_J','sky_mag_H','sky_mag_K',
                     'sky_mag_U_Moon','sky_mag_B_Moon','sky_mag_V_Moon','sky_mag_R_Moon','sky_mag_I_Moon','sky_mag_Z_Moon','sky_mag_J_Moon','sky_mag_H_Moon','sky_mag_K_Moon',
                     'skybg_airm','skybg_zod','skybg_mw','skybg_moon','zenith_ra','zenith_dec','zenith_l','zenith_b','zenith_elon','zenith_elat'])

    for i in range(nmeteo):
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

        if(i==0):
            median_pressure2m   =np.nanmedian(pressure2m[i-9:i])
            median_rhum2m       =np.nanmedian(rhum2m[i-9:i])
            median_temp2m       =np.nanmedian(temp2m[i-9:i])
            median_tempdew2m    =np.nanmedian(tempdew2m[i-9:i])
            median_winddir10m   =np.nanmedian(winddir10m[i-9:i])
            median_windspeed10m =np.nanmedian(windspeed10m[i-9:i])
            median_temp30m      =np.nanmedian(temp30m[i-9:i])
            median_winddir30m   =np.nanmedian(winddir30m[i-9:i])
            median_windspeed30m =np.nanmedian(windspeed30m[i-9:i])
            d_t            = 0
            d_pressure2m   = 0
            d_rhum2m       = 0
            d_temp2m       = 0
            d_tempdew2m    = 0
            d_winddir10m   = 0
            d_windspeed10m = 0
            d_temp30m      = 0
            d_winddir30m   = 0
            d_windspeed30m = 0
        else:
            median_pressure2m   =np.nanmedian(pressure2m[i-9:i])
            median_rhum2m       =np.nanmedian(rhum2m[i-9:i])
            median_temp2m       =np.nanmedian(temp2m[i-9:i])
            median_tempdew2m    =np.nanmedian(tempdew2m[i-9:i])
            median_winddir10m   =np.nanmedian(winddir10m[i-9:i])
            median_windspeed10m =np.nanmedian(windspeed10m[i-9:i])
            median_temp30m      =np.nanmedian(temp30m[i-9:i])
            median_winddir30m   =np.nanmedian(winddir30m[i-9:i])
            median_windspeed30m =np.nanmedian(windspeed30m[i-9:i])
            d_t            = int((jd[i]-jd[i-1])*86400)
            d_pressure2m   = pressure2m[i]  -pressure2m[i-1]
            d_rhum2m       = rhum2m[i]      -rhum2m[i-1]
            d_temp2m       = temp2m[i]      -temp2m[i-1]
            d_tempdew2m    = tempdew2m[i]   -tempdew2m[i-1]
            d_winddir10m   = winddir10m[i]  -winddir10m[i-1]
            d_windspeed10m = windspeed10m[i]-windspeed10m[i-1]
            d_temp30m      = temp30m[i]     -temp30m[i-1]
            d_winddir30m   = winddir30m[i]  -winddir30m[i-1]
            d_windspeed30m = windspeed30m[i]-windspeed30m[i-1]

        writer.writerow([
        tempo_meteo[i],d_t,jd[i],mjd[i],bjd_tdb[i],bjd_tcb[i],lst.value[i],
        year_decimal,year_frac,day_frac,day_frac_solar,year,month,dom,dow,hour_ut,hour_solar,hour_night,
        pressure2m[i],rhum2m[i],temp2m[i],tempdew2m[i],winddir10m[i],windspeed10m[i],
        temp30m[i],winddir30m[i],windspeed30m[i],
        median_pressure2m,median_rhum2m,median_temp2m,median_tempdew2m,median_winddir10m,median_windspeed10m,median_temp30m,median_winddir30m,median_windspeed30m,
        d_pressure2m,d_rhum2m,d_temp2m,d_tempdew2m,d_winddir10m,d_windspeed10m,d_temp30m,d_winddir30m,d_windspeed30m,
        sun_alt_az.alt.value[i],sun_alt_az.az.value[i],dist_sun.value[i],
        moon_alt_az.alt.value[i],moon_alt_az.az.value[i],(180-moon_phase_angle[i])/180,dist_moon.value[i],
        daylight[i],moonlight[i],solar_cycle_phase[i],
        sky_mag[i][0],sky_mag[i][1],sky_mag[i][2],sky_mag[i][3],sky_mag[i][4],sky_mag[i][5],sky_mag[i][6],sky_mag[i][7],sky_mag[i][8],
        sky_mag_wmoon[i][0],sky_mag_wmoon[i][1],sky_mag_wmoon[i][2],sky_mag_wmoon[i][3],sky_mag_wmoon[i][4],sky_mag_wmoon[i][5],sky_mag_wmoon[i][6],sky_mag_wmoon[i][7],sky_mag_wmoon[i][8],
        skybg_airm[i],skybg_zod[i],skybg_mw[i],skybg_moon[i],
        zenith_radec.ra.value[i],zenith_radec.dec.value[i],
        zenith_galac.l.value[i],zenith_galac.b.value[i],
        zenith_eclip.lon.value[i],zenith_eclip.lat.value[i]
        ])

    f_out.close()

t2=timelib.time()
print('exec time ',t2-t1,'s')
print()


# In[11]:


# determine the lenght in minutes
#tempo_tot_min=int((tempo_meteo[-1].jd-tempo_meteo[0].jd)*1440)
# create the new time vector
#newtime = tempo0 + np.linspace(0, 1, tempo_tot_min) * u.min


# In[12]:


getdimm=0
if(getdimm>0):
    query_dimm="""
SELECT start_date, ra, dec, airmass, integration, fwhm
from asm.ambient_lasilla
where start_date > '1994-01-01T00:00:00Z' and start_date< '2022-01-01T00:00:00Z'
"""
    if(FULLQUERY):
        rawframes = tapobs.search(query=query_dimm,maxrec=MAXRECORDS)
    else:
        rawframes = tapobs.search(query=query_dimm,maxrec=SMALLREC)
    tab_dimm=rawframes.to_table()
    ndimm=len(tab_dimm)
    print()
    print('----------------------')
    print("Seeing dataset")
    print('----------------------')
    print("Number of points:",ndimm)
    print()
    print("First data row")
    print(tab_dimm[0])
    print()

    #remove the last character Z from each string in the dimm dates
    t_ = [tmp_t[:-1] for tmp_t in tab_dimm['start_date']]
    tempo_dimm=Time(t_,format='isot',location=loc,scale='utc')
    # define variables from dimm
    airmass=tab_dimm['airmass']
    ra_dimm=tab_dimm['ra']
    dec_dimm=tab_dimm['dec']
    fwhm=tab_dimm['fwhm']




# In[ ]:





# In[ ]:




