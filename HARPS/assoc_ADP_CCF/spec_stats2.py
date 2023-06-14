#!/usr/bin/python3
import os
import sys
#import socket
#import psutil
#from datetime import datetime
#import copy
#import time as timecount
#from itertools import zip_longest
import pandas as pd
import numpy as np
import scipy.stats as st
from scipy.interpolate import interp1d
#import scipy.optimize
#import scipy.ndimage as ndimage
#from scipy import signal
from statsmodels.sandbox.stats import runs
from astropy.io import fits
#from astropy.timeseries import LombScargle
#from astropy import units as u
#from astropy.coordinates import Angle,SkyCoord
#from astropy.table import Table
#import healpy as hp
#import gzip
import csv
import matplotlib.pyplot as plt

def butter_lowpass(cutoff, nyq_freq, order=4):
    normal_cutoff = float(cutoff) / nyq_freq
    b, a = signal.butter(order, normal_cutoff, btype='lowpass')
    return b, a

def butter_lowpass_filter(data, cutoff_freq, nyq_freq, order=4):
    b, a = butter_lowpass(cutoff_freq, nyq_freq, order=order)
    y = signal.filtfilt(b, a, data)
    return y


def butter_highpass(cutoff, nyq_freq, order=4):
    normal_cutoff = float(cutoff) / nyq_freq
    b, a = signal.butter(order, normal_cutoff, btype='highpass')
    return b, a


def butter_highpass_filter(data, cutoff_freq, nyq_freq, order=4):
    b, a = butter_highpass(cutoff_freq, nyq_freq, order=order)
    y = signal.filtfilt(b, a, data)
    return y



tabf=pd.read_csv("adp.list")
nfiles=len(tabf)
indexf=tabf.index

#sys.exit()

wl_min_b=3783
wl_max_b=5304
wl_max_r=6912
wl_min_r=5338

with open('statistic.specs.2.NEW.csv','w') as spec_stats_file:
        specwriter=csv.writer(spec_stats_file,delimiter=',')
        specwriter.writerow(['adp','prov1','objectname','date','mjd','ra_pnt','dec_pnt','az','alt','airmass','snr','exptime','parang','specres','mag_b','mag_r','color_b_r','max_b','min_b','sum_b','mean_b','median_b','rms_b','skew_b','kurt_b','iq99_b','iq01_b','lr_m_b','lr_m_e_b','ks_v_b','ks_p_b','rt_z_b','rt_p_b','max_r','min_r','sum_r','mean_r','median_r','rms_r','skew_r','kurt_r','iq99_r','iq01_r','lr_m_r','lr_m_e_r','ks_v_r','ks_p_r','rt_z_r','rt_p_r'
        'rms_bf','skew_bf','kurt_bf','rms_rf','skew_rf','kurt_rf'])
        spec_stats_file.close()


for i in range(nfiles):
        stat_spec=[]
        fname=tabf.filename[i]
        if(i>0):
                fsize=os.path.getsize(fname)
                if(fsize>0):
                    print( 'INFO ',i,',',fname,', file size =', fsize)
                    hdu=fits.open(fname,ignore_missing_end=True)
                    hdr0    = hdu[0].header
                    hdr1    = hdu[1].header
                    spec    = hdu[1].data
                    wl=0.0
                    flux=0.0
                    wl      = spec.field('WAVE')
                    flux    = spec.field('FLUX')
                    npoints = flux.size
                    adp=fname
                    ra_pnt=hdr0["RA"]
                    dec_pnt=hdr0["DEC"]
                    objectname=hdr0["OBJECT"]
                    snr=hdr0["SNR"]
                    exptime=hdr0["EXPTIME"]
                    mjd=hdr0["MJD-OBS"]
                    date=hdr0["DATE-OBS"]
                    airmass=hdr0["HIERARCH ESO TEL AIRM START"]
                    specres=hdr0["SPEC_RES"]
                    prov1=hdr0["PROV1"]
                    parang=hdr0["HIERARCH ESO TEL PARANG START"]
                    az=hdr0["HIERARCH ESO TEL AZ"]
                    alt=hdr0["HIERARCH ESO TEL ALT"]
                    az=np.mod(az,360.0)

                    if(npoints>150000):
                            #print(npoints)
                            j1=0
                            j2=0
                            j3=0
                            j4=0

                            for j in range(0,300):
                                    if(wl[0][j-1]<wl_min_b and wl[0][j]>=wl_min_b):
                                            j1=j

                            for j in range(152100,152400):
                                    if(wl[0][j-1]<wl_max_b and wl[0][j]>=wl_max_b):
                                            j2=j

                            for j in range(155400,155700):
                                    if(wl[0][j-1]<wl_min_r and wl[0][j]>=wl_min_r):
                                            j3=j

                            for j in range(npoints-300,npoints-1):
                                    if(wl[0][j-1]<wl_max_r and wl[0][j]>=wl_max_r):
                                            j4=j

                            if(j1>0 and j2>j1):
                                x=flux[0][j1:j2]
                                wl_b=wl[0][j1:j2]
                                max_b       = np.nanmax(x)
                                min_b       = np.nanmin(x)
                                sum_b       = np.nansum(x)
                                mean_b=np.nanmean(x)
                                median_b=np.nanmedian(x)
                                rms_b=np.nanstd(x,ddof=1)
                                skew_b=st.skew(x,nan_policy='omit')
                                kurt_b=st.kurtosis(x,nan_policy='omit')
                                iq99_b=st.iqr(x,rng=[50,99],nan_policy='omit')
                                iq01_b=st.iqr(x,rng=[1,50],nan_policy='omit')
                                mag_b=-2.5*np.log10(sum_b)
                                res = st.linregress(wl_b,x)
                                lr_m_b = res.slope
                                lr_m_e_b = res.stderr
                                ks_v_b,ks_p_b=st.kstest(x,'norm')
                                rt_z_b,rt_p_b=runs.runstest_1samp(x,cutoff='median',correction=False)#Runstest


                            coarse_resolution = 100.0  # Angstrom
                            nnn=int(round((wl_b[-1]-wl_b[0])/coarse_resolution))
                            nwl_b = np.linspace(wl_b[0], wl_b[-1], nnn)
                            spectrum_b_LR = np.zeros_like(nwl_b)
                            for i in range(nnn):
                                    start_idx = np.searchsorted(wl_b, nwl_b[i] - coarse_resolution / 2)
                                    end_idx = np.searchsorted(wl_b, nwl_b[i] + coarse_resolution / 2)
                                    spectrum_b_LR[i] = np.sum(x[start_idx:end_idx])/(end_idx-start_idx)

                            f = interp1d(nwl_b, spectrum_b_LR, kind='linear')
                            spectrum_b_envelope = f(wl_b)
                            spectrum_b_flat_0 = x / spectrum_b_envelope
                            spec_median_b = np.nanmedian(spectrum_b_flat_0)
                            spectrum_b_flat = spectrum_b_flat_0/spec_median_b

                            rms_bf=np.nanstd(spectrum_b_flat,ddof=1)
                            skew_bf=st.skew(spectrum_b_flat,nan_policy='omit')
                            kurt_bf=st.kurtosis(spectrum_b_flat,nan_policy='omit')


                            if(j3>j2 and j3<j4):
                                x=flux[0][j3:j4]
                                wl_r=wl[0][j3:j4]
                                max_r       = np.nanmax(x)
                                min_r       = np.nanmin(x)
                                sum_r       = np.nansum(x)
                                mean_r=np.nanmean(x)
                                median_r=np.nanmedian(x)
                                rms_r=np.nanstd(x,ddof=1)
                                skew_r=st.skew(x,nan_policy='omit')
                                kurt_r=st.kurtosis(x,nan_policy='omit')
                                iq99_r=st.iqr(x,rng=[50,99],nan_policy='omit')
                                iq01_r=st.iqr(x,rng=[1,50],nan_policy='omit')
                                mag_r=-2.5*np.log10(sum_r)
                                res = st.linregress(wl_r,x)
                                lr_m_r = res.slope
                                lr_m_e_r = res.stderr
                                ks_v_r,ks_p_r=st.kstest(x,'norm')
                                rt_z_r,rt_p_r=runs.runstest_1samp(x,cutoff='median',correction=False)#Runstest

                            color_b_r=mag_b-mag_r

                            nnn=int(round((wl_r[-1]-wl_r[0])/coarse_resolution))
                            nwl_r = np.linspace(wl_r[0], wl_r[-1], nnn)
                            spectrum_r_LR = np.zeros_like(nwl_r)
                            for i in range(len(nwl_r)):
                                    start_idx = np.searchsorted(wl_r, nwl_r[i] - coarse_resolution / 2)
                                    end_idx = np.searchsorted(wl_r, nwl_r[i] + coarse_resolution / 2)
                                    spectrum_r_LR[i] = np.sum(x[start_idx:end_idx])/(end_idx-start_idx)
                            f = interp1d(nwl_r, spectrum_r_LR, kind='linear')
                            spectrum_r_envelope = f(wl_r)
                            spectrum_r_flat_0 = x / spectrum_r_envelope
                            spec_median_r = np.nanmedian(spectrum_r_flat_0)
                            spectrum_r_flat = spectrum_r_flat_0/spec_median_r

                            rms_rf=np.nanstd(spectrum_r_flat,ddof=1)
                            skew_rf=st.skew(spectrum_r_flat,nan_policy='omit')
                            kurt_rf=st.kurtosis(spectrum_r_flat,nan_policy='omit')


                            #plt.plot(wl_b, spectrum_b_flat)
                            #plt.plot(wl_r, spectrum_r_flat)
                            #plt.plot(nwl_b, spectrum_b_LR)
                            #plt.plot(nwl_r, spectrum_r_LR)
                            #plt.xlabel('Wavelength (Angstrom)')
                            #plt.ylabel('Flux')
                            #plt.show()

                            if(j1>0 and j2>j1 and j3>j2 and j3<j4):
                                stat_spec=[adp,prov1,objectname,date,mjd,ra_pnt,dec_pnt,az,alt,airmass,snr,exptime,parang,specres,mag_b,mag_r,color_b_r,max_b,min_b,sum_b,mean_b,median_b,rms_b,skew_b,kurt_b,iq99_b,iq01_b,lr_m_b,lr_m_e_b,ks_v_b,ks_p_b,rt_z_b,rt_p_b,max_r,min_r,sum_r,mean_r,median_r,rms_r,skew_r,kurt_r,iq99_r,iq01_r,lr_m_r,lr_m_e_r,ks_v_r,ks_p_r,rt_z_r,rt_p_r,
                                rms_bf,skew_bf,kurt_bf,rms_rf,skew_rf,kurt_rf]
                                with open('statistic.specs.2.NEW.csv','a') as spec_stats_file:
                                    specwriter=csv.writer(spec_stats_file,delimiter=',')
                                    specwriter.writerow(stat_spec)
                                    spec_stats_file.close()
                            else:
                                print( 'ERR 3, ',fname,', j1 =', j1)
                                print( 'ERR 3, ',fname,', j2 =', j2)
                                print( 'ERR 3, ',fname,', j3 =', j3)
                                print( 'ERR 3, ',fname,', j4 =', j4)

                            #sys.exit()

                    else:
                        print( 'ERR 2, ',fname,', num. points =', fsize)
                else:
                        print( 'ERR 1, ',fname,', file size =', fsize)

#        stat_spec.append([adp,prov1,objectname,date,mjd,ra_pnt,dec_pnt,az,alt,airmass,snr,exptime,parang,specres,mag_b,mag_r,color_b_r,max_b,min_b,sum_b,mean_b,median_b,rms_b,skew_b,kurt_b,iq99_b,iq01_b,lr_m_b,lr_m_e_b,ks_v_b,ks_p_b,rt_z_b,rt_p_b,max_r,min_r,sum_r,mean_r,median_r,rms_r,skew_r,kurt_r,iq99_r,iq01_r,lr_m_r,lr_m_e_r,ks_v_r,ks_p_r,rt_z_r,rt_p_r])




#with open('stats.csv','w') as spec_stats_file:
#        specwriter=csv.writer(spec_stats_file,delimiter=',')
#        specwriter.writerow(['adp','prov1','objectname','date','mjd','ra_pnt','dec_pnt','az','alt','airmass','snr','exptime','parang','specres','mag_b','mag_r','color_b_r','max_b','min_b','sum_b','mean_b','median_b','rms_b','skew_b','kurt_b','iq99_b','iq01_b','lr_m_b','lr_m_e_b','ks_v_b','ks_p_b','rt_z_b','rt_p_b','max_r','min_r','sum_r','mean_r','median_r','rms_r','skew_r','kurt_r','iq99_r','iq01_r','lr_m_r','lr_m_e_r','ks_v_r','ks_p_r','rt_z_r','rt_p_r'])
#        for i in range(len(stat_spec)):
#                specwriter.writerow(stat_spec[i])

