#!/usr/bin/python3
# versione con modifiche suggerite dal referee

#import matplotlib.pyplot as plt

import os
import sys
import numpy as np
import gzip
import csv
import pandas as pd
import scipy.stats as st
import scipy.optimize
import scipy.ndimage as ndimage
from scipy import interpolate
from scipy import signal
from scipy.stats import sigmaclip
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import angular_separation
import copy
import matplotlib.pyplot as plt


d2r = np.pi/180.0
r2d = 180.0/np.pi


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



def func_clip(x,method,Ltrh,Utrh,wlen):
# fa un sigma clipping con sostituzione
# prende la serie originale e trova i valori oltre soglie di rms
# questi valori vengono sostituiti dalla media aritmetica del punto precedente e seguente
# se i punti sono all'inizio o fine delle sarie si calcola la media sui 2 punti seguenti o precedenti
# x_ e' la serie di dati
# method : calcola i limiti superiori o inferiori con la media (1) mediana (2)
# Ltrh : limite inferiore in sigma
# Utrh : limite superiore in sigma

# calcola i valori utili
        x_rms       = np.nanstd(x,ddof=1)
        if(method==1):
                 x_ave = np.nanmean(x)
        elif(method==2):
                 x_ave = np.nanmedian(x)
        Lbound = x_ave - x_rms*Ltrh
        Ubound = x_ave + x_rms*Utrh
        nout=0
        nin=0
        xnew=[]
        n=len(x)
        mask_clip=np.zeros(n)
# fa il loop su tutti i valori della serie e li sostituisce
        for i in range(n):
                imin=max(0,i-wlen)
                imax=min(i+wlen,n)
                Lmean   = np.nanmean(x[imin:imax])
                Lmedian = np.nanmedian(x[imin:imax])
                Lstd    = np.nanstd(x[imin:imax])
                if(method==1):
                        x1      = (x[i]-Lmean)/Lstd
                elif(method==2):
                        x1      = (x[i]-Lmedian)/Lstd
                if( x1>=Lbound and x1<=Ubound):
                        nin=nin+1
                        xnew.append(x[i])
                        mask_clip[i]=0
                else:
                        nout=nout+1
                        mask_clip[i]=1

        mean   = np.nanmean(xnew)
        median = np.nanmedian(xnew)
        std    = np.nanstd(xnew)
        return(xnew,mask_clip,nin,nout,mean,median,std)





def moving_mode(x,wlen):
        xnew=[]
        n=len(x)
        for i in range(n):
                imin=max(0,i-wlen)
                imax=min(i+wlen,n)
                h = np.histogram(x[imin:imax],bins=31)
                k=h[0][:].argmax() # index of the max
                localmode=h[1][k]  # value of the peak
                xnew.append(localmode)
        return(xnew)



def read_spec(filename):
    '''Read a HARPSN spectrum

    Parameters
    ----------
    filename : string
    name of the fits file with the data

    Returns
    -------
    wavelength : np.ndarray
    wavelength (in Ang)
    flux : np.ndarray
    flux (in erg/s/cm**2)
    date_obs : string
    time of observation
    '''
    sp = fits.open(filename)
    header = sp[0].header

    wcs = WCS(header)
    #make index array
    index = np.arange(header['NAXIS1'])

    wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
    wavelength = wavelength.flatten()
    flux = sp[0].data

    ra             = header['RA-DEG']
    dec            = header['DEC-DEG']
    airmass        = header['AIRMASS']
    elevation      = header['EL']*180.0/np.pi
    azimuth        = header['HIERARCH TNG VTRK TRK AZENC']*180.0/np.pi
    date           = header['Date-OBS']
    mjd            = header['MJD-OBS']
    exptime        = header['EXPTIME']
    berv           = header['HIERARCH TNG DRS BERV']
    fwhmx          = header['HIERARCH TNG AG ACQ FWHMX']
    fwhmy          = header['HIERARCH TNG AG ACQ FWHMY']
    posx           = header['HIERARCH TNG AG ACQ POSX']
    posy           = header['HIERARCH TNG AG ACQ POSY']
    name           = header['HIERARCH TNG OBS TARG NAME']

    return (wavelength,flux,name,date,mjd,ra,dec,airmass,elevation,azimuth,exptime,berv,fwhmx,fwhmy,posx,posy)



def read_ccf(filename):
    '''Read a HARPSN spectrum CCF

    Parameters
    ----------
    filename : string
    name of the fits file with the data

    Returns
    -------
    wavelength : np.ndarray
    wavelength (in Ang)
    flux : np.ndarray
    flux (in erg/s/cm**2)
    date_obs : string
    time of observation
    '''
    ccf = fits.open(filename)
    header = ccf[0].header

    wcs = WCS(header)

    rv             = header['HIERARCH TNG DRS CCF RVC']
    rv_e           = header['HIERARCH TNG DRS CCF NOISE']

    return (rv,rv_e)


def norm_spec(x,y):

        order=3
        pol,residuals,rank,svd,rcond = np.polyfit(x,y,order,full=True,cov=True)
        y_smooth = np.polyval(pol,x)
        y_flat      = y / y_smooth

        #fig, subfig = plt.subplots(1)
        ####subfig.plot(x,y)
        ####subfig.plot(x,y_smooth)
        #subfig.plot(x,y_flat)
        #plt.show()


        return(y_flat)










#w1=3875
#w2=4500

#w1=4500
#w2=5500

w1=5500
w2=6450

#w1=6450
#w2=6908

#w1=3875
#w2=6908

#w1=5400
#w2=6450

#w1=3900
#w2=6900


qe_data=pd.read_csv("harpn_qe.csv")
finterp_qe=interpolate.interp1d(qe_data["lambda"],qe_data["qe"],fill_value="extrapolate")

extco_data=pd.read_csv("ext_coeff_n1.csv")
finterp_extcoeff=interpolate.interp1d(extco_data["lambda"],extco_data["extcoeff"],fill_value="extrapolate")



s95=pd.read_csv("spencer1995.csv")
wspen=s95["lambda"]
fspen=s95["ganymede"]





lfiles=pd.read_csv("lista_ref_n1.csv")
nfiles=len(lfiles)



######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
# telluric star
i=0
fs1d=lfiles['s1d'][i]
fccf=lfiles['ccf'][i]
wavelength,flux,name,date,mjd,ra,dec,airmass,elevation,azimuth,exptime,berv,fwhmx,fwhmy,posx,posy=read_spec(fs1d)
print(lfiles['s1d'][i],name,date)
print(elevation)
el0=elevation
az0=azimuth
mjd0=mjd
am0=airmass
#rv0,rv_e0=read_ccf(fccf) ## wrong mask
#print(rv0)
rv0=-11.0
rv_e0=4.0
vrad0=rv0-berv
wl=wavelength[(wavelength>w1) & (wavelength<w2)]
ff      =flux[(wavelength>w1) & (wavelength<w2)]

dw=300
flusso=[]
lungon=[]
j=-1
for w in range(w1,w2,dw):
        j=j+1
        fl_=np.sum(ff[(wl>=w) & (wl<w+dw)])
        lungon.append(w+dw/2)
        flusso.append(fl_)
finterp_flusso=interpolate.interp1d(lungon,flusso,fill_value="extrapolate")
flusso_interp=finterp_flusso(wl)
flusso_interp=flusso_interp/np.max(flusso_interp)



spec_tellstar=norm_spec(wl,ff)
wl_tell = wl

#spec_tellstar=ff/flusso_interp

#fig, subfig = plt.subplots(1)
#subfig.plot(wl,ff_/5e4)
#subfig.plot(wl,diomerda)
#subfig.plot(wl,qe_interp)
#subfig.plot(wl,extcoeff_interp)
#plt.show()
#sys.exit()






######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
# solar analog
i=1
fs1d=lfiles['s1d'][i]
fccf=lfiles['ccf'][i]
wavelength,flux,name,date,mjd,ra,dec,airmass,elevation,azimuth,exptime,berv,fwhmx,fwhmy,posx,posy=read_spec(fs1d)
print(lfiles['s1d'][i],name,date)
print(elevation)
el1=elevation
az1=azimuth
mjd1=mjd
am1=airmass
rv1,rv_e1=read_ccf(fccf)
rv1=-5.724
rv_e1=0.0003
wl=wavelength[(wavelength>w1) & (wavelength<w2)]
ff      =flux[(wavelength>w1) & (wavelength<w2)]


dw=300
flusso=[]
lungon=[]
j=-1
for w in range(w1,w2,dw):
        j=j+1
        fl_=np.sum(ff[(wl>=w) & (wl<w+dw)])
        lungon.append(w+dw/2)
        flusso.append(fl_)
finterp_flusso=interpolate.interp1d(lungon,flusso,fill_value="extrapolate")
flusso_interp=finterp_flusso(wl)
flusso_interp=flusso_interp/np.max(flusso_interp)

#k=ff[:].argmax() # index of the max
#localmax=wl[k]  # value of the peak
#print(k,wl[k],ff[k])
#k=flusso_interp[:].argmax() # index of the max
#localmax=wl[k]  # value of the peak
#print(k,wl[k],ff[k])
#sys.exit()



#print('qui')
#splina = interpolate.BSpline(wl,ff,s=10)
#print('quo')

#fig, subfig = plt.subplots(1)
#subfig.plot(wl,ff)
#subfig.plot(wl,ff/flusso_interp)
#plt.show()






spec_solanalog=norm_spec(wl,ff)
#print(np.mean(spec_tellstar),np.mean(spec_solanalog))

#####spec_solanalog=ff/flusso_interp
#print(np.mean(spec_tellstar),np.mean(spec_solanalog))

vrad1=rv1-berv

wl_solar=wl

#sys.exit()

calcola_rv=0
if(calcola_rv==1):
        vr=[]
        s=[]
        rv_min=28.6-5
        rv_span=10
        nv=20001
        rv_step=rv_span/nv
        l1=6270
        l2=6290
        wl_tell1      =wl_tell[(wl_tell>l1) & (wl_tell<l2)]
        spec_tellstar1=spec_tellstar[(wl_tell>l1) & (wl_tell<l2)]
        wl_solar1     =wl_solar[(wl_solar>l1) & (wl_solar<l2)]
        spec_solar1   =spec_solanalog[(wl_solar>l1) & (wl_solar<l2)]


        for j in range(nv):
                velrad=rv_min+j*rv_step
                vr.append(velrad)
                wl_o2=wl_solar1*(1.0+velrad/299792.458)
                finterp=interpolate.interp1d(wl_tell1,spec_tellstar1,fill_value="extrapolate")
                spec_tellstar_interp=finterp(wl_o2)
                spec_solatell=spec_solar1-spec_tellstar_interp
                s.append(np.std(spec_solatell))
                print(j,s[j])
        k=s.index(np.min(s)) # indice della RV che minimizza i residui
        print(vr[k])
        # miglior valore 28.6175657
        fig, subfig = plt.subplots(1)
        subfig.plot(vr,s)
        plt.show()


#print(np.mean(spec_tellstar),np.mean(spec_solanalog))

#sys.exit()

rv_shift=28.692245387730615
#rv_shift=vr[k]

wl_o2=wl_solar*(1.0+rv_shift/299792.458)

finterp=interpolate.interp1d(wl_tell,spec_tellstar,fill_value="extrapolate")
spec_tellstar_interp=finterp(wl_o2)

spec_solatell=spec_solanalog-spec_tellstar_interp+1.0

l1=6275
l2=6280
#plt.figure(figsize=(18,12))
#plt.xlim(l1,l2)
#plt.plot(wl_solar,spec_solanalog)
##plt.plot(wl_tell,spec_tellstar)
#plt.plot(wl_o2,spec_tellstar_interp)
#plt.plot(wl_solar,spec_solatell)
#plt.show()
#sys.exit()


#spec_solatell=spec_solanalog-spec_tellstar_interp+1.0



#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################

#f_1=open('spec_ref.csv','wt',newline='')
#scrivi = csv.writer(f_1,delimiter=',')
#scrivi.writerow(['wl_obs','wl_o2','spec_solatell','spec_tell','spec_solanalog'])
#for i in range(len(wl)):
#        scrivi.writerow([wl[i],wl_o2[i],spec_solatell[i],spec_tellstar_interp[i],spec_solanalog[i]])
#f_1.close()



logobs=pd.read_csv("observation_n1+ephem.log.csv")
nfiles=len(logobs)
f_gplt=open('figure_spec.plt','a',newline='')

avg_spec=np.zeros(len(spec_solatell))


nfiles=121
#nfiles=1

for i in range(nfiles):
        fs1d=logobs['fs1d'][i]
        fccf=logobs['fccf'][i]
        fout=logobs['fs1d'][i][0:29]+'.spec.csv'
        wavelength,flux_,name,date,mjd,ra,dec,airmass,elevation,azimuth,exptime,berv,fwhmx,fwhmy,posx,posy=read_spec(fs1d)
        print(i,fs1d,name,date,elevation)
        rv,rv_e=read_ccf(fccf)
        angsep=angular_separation(az0*d2r,el0*d2r,azimuth*d2r,elevation*d2r)*r2d
        dt=(mjd-mjd0)*24
        if(i==0):
                rv0gan=rv
        drv=rv-rv0gan
        qe_interp=finterp_qe(wavelength)
        flux=flux_/qe_interp
        wl=wavelength[(wavelength>w1) & (wavelength<w2)]
        ff      =flux[(wavelength>w1) & (wavelength<w2)]
        #rdot=logobs['rdot'][i]
        #deldot=logobs['deldot'][i]
        #vrad=rv-berv#+rdot+deldot
        #sunshift=rv-berv-rdot-deldot


        dw=300
        flusso=[]
        lungon=[]
        j=-1
        for w in range(w1,w2,dw):
                j=j+1
                fl_=np.sum(ff[(wl>=w) & (wl<w+dw)])
                lungon.append(w+dw/2)
                flusso.append(fl_)
        finterp_flusso=interpolate.interp1d(lungon,flusso,fill_value="extrapolate")
        flusso_interp=finterp_flusso(wl)
        flusso_interp=flusso_interp/np.max(flusso_interp)

        # spec_gan0=ff/flusso_interp ### WRONG

        spec_gan0 =ff

        spec_gan1=norm_spec(wl,ff)

        #plt.plot(wl,ff)
        #plt.plot(wl,spec_gan0)
        #plt.plot(wl,spec_gan1)
        #plt.show()
        #sys.exit()


        #spec_gan1=spec_gan0-spec_tellstar+1.0 # sottraggo la tellurica
        #spec_gan1=spec_gan0/spec_tellstar # divido per la tellurica
        if(i==0):
                spec_gan1_ref=spec_gan1


        calcola_rv=0 # tra spettri consecutivi di Ganimede
        if(calcola_rv==1):
                vr=[]
                s=[]
                nv=401
                rv_step=0.01
                rv_guess=-0.5
                for j in range(nv):
                        velrad=rv_guess+j*rv_step
                        vr.append(velrad)
                        wl_gh=wl*(1.0+velrad/299792.458)
                        finterpa=interpolate.interp1d(wl,spec_gan0,fill_value="extrapolate")
                        spec_gan1_interp=finterpa(wl_gh)
                        spec_g1=spec_gan1_ref-spec_gan1_interp
                        s.append(np.std(spec_g1))
                k=s.index(np.min(s)) # indice della RV che minimizza i residui
                print(i,fout,vr[k],rv-rv0gan,vr[k]-(rv-rv0))
                # miglior valore 28.6175657
                #fig, subfig = plt.subplots(1)
                #subfig.plot(vr,s)
                #plt.show()


        #wl0_gan=wl*(1.0+(rv1-rv)/299792.458) # creo le lunghezze d'onda relative alla differenza tra HD221627 e Ganimede
        #finterp1=interpolate.interp1d(wl,spec_solatell,fill_value="extrapolate")
        #spec_solatell_interp=finterp1(wl0_gan) # interpolo lo spettro di HD221627 per riportarlo alla rv di Ganimede
        #spec_gan2=spec_gan1-spec_solatell_interp+1.0



        calcola_rv=1
        if(calcola_rv==1):
                vr=[]
                s=[]
                nv=201
                rv_min=-25
                rv_span=10
                rv_step=rv_span/nv
                for j in range(nv):
                        velrad=rv_min+j*rv_step
                        vr.append(velrad)
                        wl_gh=wl*(1.0+velrad/299792.458)
                        finterp1=interpolate.interp1d(wl,spec_solatell,fill_value="extrapolate")
                        spec_solatell_interp=finterp1(wl_gh)
                        spec_g2=spec_gan1-spec_solatell_interp
                        ##print(np.mean(spec_solatell_interp),np.mean(spec_gan1))
                        s.append(np.std(spec_g2))
                        l1=6275
                        l2=6285
                        #plt.plot(wl_gh,spec_gan1)
                        #plt.plot(wl_gh,spec_solatell_interp)
                        #plt.xlim(l1,l2)
                        #plt.show()
                        #plt.close
                k=s.index(np.min(s)) # indice della RV che minimizza i residui
                print(i,fout,vr[k],rv1-rv,vr[k]-(rv1-rv))
                #fig, subfig = plt.subplots(1)
                #subfig.plot(vr,s)
                #plt.show()
        #sys.exit()

        wl0_gan=wl*(1.0+vr[k]/299792.458)
        finterp1=interpolate.interp1d(wl,spec_solatell,fill_value="extrapolate")
        spec_solatell_interp=finterp1(wl0_gan)
        spec_gan2=spec_gan1-spec_solatell_interp+1.0


        # rimuove la tellurica
        calcola_rv=1
        if(calcola_rv==1):
                vr=[]
                s=[]
                rv_min=-1
                rv_span=2
                nv=201
                rv_step=rv_span/nv
                l1=6270
                l2=6290
                wl_tell1      =wl_tell[(wl_tell>l1) & (wl_tell<l2)]
                spec_tellstar1=spec_tellstar[(wl_tell>l1) & (wl_tell<l2)]
                wl_gan22     =wl[(wl>l1) & (wl<l2)]
                spec_gan22   =spec_gan2[(wl>l1) & (wl<l2)]


                for j in range(nv):
                        velrad=rv_min+j*rv_step
                        vr.append(velrad)
                        wl_o2=wl_gan22*(1.0+velrad/299792.458)
                        finterp=interpolate.interp1d(wl_tell1,spec_tellstar1,fill_value="extrapolate")
                        spec_tellstar_interp=finterp(wl_o2)
                        spec_gantell=spec_gan22-spec_tellstar_interp
                        s.append(np.std(spec_gantell))
                        #print(j,s[j])
                k=s.index(np.min(s)) # indice della RV che minimizza i residui
                print(vr[k])
                #fig, subfig = plt.subplots(1)
                #subfig.plot(vr,s)
                #plt.show()


        #sys.exit()

        rv_shift=vr[k]

        wl_o2=wl*(1.0+rv_shift/299792.458)

        finterp=interpolate.interp1d(wl_tell,spec_tellstar,fill_value="extrapolate")
        spec_tellstar_interp=finterp(wl_o2)

        spec_gantell=spec_gan2-spec_tellstar_interp+1.0
        #plt.plot(wl,spec_gantell,'.')
        #plt.show()

        #sys.exit()







        wl0_gan_air=wl*(1.0+(-rv)/299792.458) # creo le lunghezze d'onda relative alla rv di Ganimede
        finterpg=interpolate.interp1d(wl,spec_gantell,fill_value="extrapolate")
        spec_gan1_shifted=finterpg(wl0_gan_air) # interpolo lo spettro di Ganimede per riportarlo al rv=0

        #import matplotlib.pyplot as plt
        #fig, subfig = plt.subplots(nrows=2,ncols=1,figsize=(18,12))
        #subfig[0].plot(wl,spec_gan1,'-')
        #subfig[0].plot(wl,spec_solatell_interp,'-')
        #subfig[0].set_xlim([w1,w2])
        #subfig[1].plot(wl,spec_gan1-spec_solatell_interp,'-')
        #subfig[1].set_xlim([w1,w2])
        #plt.show()

        #sys.exit()


        spec_gan3=ndimage.uniform_filter1d(spec_gantell,400,mode='mirror') # moving average
        #spec_gan3=ndimage.gaussian_filter1d(spec_gan2,500,mode='mirror') # gaussian kernel

        if(i==0):
                spec_gan2_ref=spec_gantell

        spec_gan2_evol=spec_gantell-spec_gan2_ref

        scrivi_file=1
        if(scrivi_file==1):
                f_sp=open(fout,'wt',newline='')
                scrivi = csv.writer(f_sp,delimiter=',')
                scrivi.writerow(['wl0','wl','spec_gan','spec_gan_norm','spec_gan_tell','spec_gan_solatell','spec_gan_solatell_smooth','spec_gan_evol','spec_tell','spec_solatell_interp','spec_solatell','spec_solanalog'])
                for j in range(len(ff)):
                        scrivi.writerow([wl0_gan[j],wl[j],ff[j],spec_gan0[j],spec_gan1[j],spec_gantell[j],spec_gan3[j],spec_gan2_evol[j],spec_tellstar[j],spec_solatell_interp[j],spec_solatell[j],spec_solanalog[j]])
                f_sp.close()

        #avg_spec=avg_spec+spec_gan1
        avg_spec=avg_spec+spec_gantell

        #vs=0.2*i
        #print('"'+fout+'" u 2:($6+',vs,') w l,\\')

#        print(i,fout)


print('*******************************')
#sys.exit()

avg_spec=avg_spec/nfiles
refspec = avg_spec
nspec=len(avg_spec)
refspec_movingmode = moving_mode(refspec,401)
finterp=interpolate.interp1d(wspen,fspen,fill_value="extrapolate")
fspencer95=finterp(wl)

#for window in range(11,110,10):
#        refspec_movingmode = moving_mode(refspec,window)
#        mean=np.mean(refspec_movingmode)
#        rms=np.std(refspec_movingmode)
#        chi2=np.sum((refspec_movingmode/mean-fspencer95)**2)/nspec
#        print(window,chi2,rms,mean)
#sys.exit()

for window in range(500,3001,500):
        ofile='spec_ref_wlen_'+str(window)+'.csv'
        print('WINDOW LENGTH',window,ofile)
        #

        refspec_ma1=ndimage.uniform_filter1d(refspec_movingmode,window,mode='mirror') # moving average
        refspec_ma2=ndimage.gaussian_filter1d(refspec_movingmode,window,mode='mirror') # gaussian kernel
        resolution=0.01 # AA per pixel
        sample_rate = 1/resolution
        cutoff_frequency = 1/(window*resolution)
        refspec_ma3=butter_lowpass_filter(refspec_movingmode, cutoff_frequency, sample_rate/2)


        #refspec_clipped, refspec_mask_clipped, refspec_clipped_nclipped, refspec_clipped_ngood, refspec_clipped_mean, refspec_clipped_median, refspec_clipped_std = func_clip(refspec,2,2.5,1.5,2000)
        #print(refspec_clipped_nclipped)
        #wl_new = wl[refspec_mask_clipped==0]
        #finterp2=interpolate.interp1d(wl_new,refspec_clipped,fill_value="extrapolate")
        #refspec_interp=finterp2(wl)
        #refspec_ma1=ndimage.uniform_filter1d(refspec_interp,1200,mode='mirror') # moving average
        #refspec_ma2=ndimage.gaussian_filter1d(refspec_interp,1000,mode='mirror') # gaussian kernel
        #resolution=0.01 # AA per pixel
        #sample_rate = 1/resolution
        #cutoff_frequency = 1/(1000*resolution) #  8 points
        #refspec_ma3=butter_lowpass_filter(refspec_interp, cutoff_frequency, sample_rate/2)



        # ma1 = moving avg
        # ma2 = gaussian kernel
        # ma3 = low pass

        #splint=interpolate.splrep(wl,refspec,s=0)
        #refspec_spline=interpolate.splev(wl, splint, der=0)
        f_rsp=open(ofile,'wt',newline='')
        scrivi = csv.writer(f_rsp,delimiter=',')
        scrivi.writerow(['lambda','refspec','refspec_movingmode','refspec_ma','refspec_gk','refspec_lp','Ganymede','HD221627','Spencer95'])
        for i in range(nspec):
                scrivi.writerow([wl[i],refspec[i],refspec_movingmode[i],refspec_ma1[i],refspec_ma2[i],refspec_ma3[i],avg_spec[i],spec_solatell_interp[i],fspencer95[i]])
        f_rsp.close()

        plt.figure(figsize=(12,9))
        plt.xlim(wl[0],wl[nspec-1])
        plt.ylim(0.97,1.01)
        ##plt.plot(wl,refspec)
        #plt.plot(wl,refspec_movingmode) # blue
        plt.plot(wl,refspec_ma1)        # green
        #plt.plot(wl,refspec_ma2)        # red 
        #plt.plot(wl,refspec_ma3)        # orange
        plt.plot(wl,fspencer95-0.005)   # magenta
        plt.savefig('ref_spec_wlen_'+str(window)+'.png', format='png')
        #plt.show()
        plt.close()




