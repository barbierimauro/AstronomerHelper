#!/scratch/mbarbier/miniconda3/bin/python3
#
# (c) mauro barbieri 2023
# maurobarbieri.science@gmail.com
#

#********************************************************************************
#
#  PROGRAM SETUP
#
#********************************************************************************

import os
import sys
import glob
import warnings
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal, interpolate, stats
from scipy.signal import correlate, find_peaks, peak_widths
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import astropy.units as u
from astropy.io import fits
from astropy import constants as const
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_body, get_sun, get_moon
from astropy.time import Time


# Filter the warning
logging.basicConfig(filename='error.log', level=logging.ERROR)
warnings.filterwarnings('ignore', category=UserWarning, append=True)


#********************************************************************************
#
#  FUNCTIONS
#
#********************************************************************************

#placeholder for new functions

# ------------------------------  BEGIN FUNCTION ------------------------------ #
#def funcname():
# ------------------------------  END   FUNCTION ------------------------------ #





# ------------------------------  BEGIN FUNCTION ------------------------------ #
def fit_absorption_line2(wavelengths, intensities, wc,ws):
    line_center0 = wc
    wl_span = ws
    mask = (wavelengths > (line_center0 - wl_span)) & (wavelengths < (line_center0 + wl_span))
    wl=wavelengths[mask]
    mf=np.mean(intensities[mask])
    flux_i = 1e6 - intensities[mask]
    flux = intensities[mask]
    npt=len(wavelengths[mask])
    
    
    smoothed_flux_i = gaussian_filter(flux_i, sigma=5)
    smoothed_flux0_i = gaussian_filter(flux_i, sigma=100)
    residuals_i = flux_i - smoothed_flux_i
    rms_i = np.sqrt(np.mean(residuals_i**2))
    threshold_i = 5 * rms_i
    peaks_i, properties_i = find_peaks(smoothed_flux_i, width=11, height=0.01)
    filtered_peaks_i = []
    for i, peak_i in enumerate(peaks_i):
        width_i = properties_i['widths'][i]
        peak_value_i = properties_i['peak_heights'][i]
        local_background_i = smoothed_flux0_i[peak_i]
        local_rms_i = rms_i  
        if peak_value_i > local_background_i + 3 * local_rms_i and width_i > 11:
            filtered_peaks_i.append((peak_i, width_i, peak_value_i))
    
    # Sort filtered_peaks by width in descending order
    filtered_peaks_i.sort(key=lambda x: x[1], reverse=True)    
    for i, (peak_i, width_i, peak_value_i) in enumerate(filtered_peaks_i, start=1):
        print(f"{wl[peak_i]:7.2f},{width_i * 0.01:7.3f},{1-(1e6 - peak_value_i)/mf:5.2f}")

    peak_indices_i = [peak_info[0] for peak_info in filtered_peaks_i]
    
    print('n peaks', len(filtered_peaks_i))

#    if filtered_peaks:        
#        const_guess = np.mean(flux)
#        amplitude_guess = filtered_peaks[0][2]
#        mean_guess = wl[filtered_peaks[0][0]]
#        stddev_guess = filtered_peaks[0][1]
#    else:        
    const_guess = np.mean(flux)
    amplitude_guess = np.min(flux) - const_guess
    mean_guess = line_center0
    stddev_guess = 0.2

    #fit Gaussian
    ERRORE=0
    try:
        popt_gaussian, pcov = curve_fit(negative_gaussian, wl, flux, p0=[const_guess, amplitude_guess, mean_guess, stddev_guess], maxfev=2000)
    except RuntimeError as e:
        ERORRE = 1
        popt_gaussian= [np.nan for _ in range(4)]
        pcov = [[np.nan for _ in range(4)] for _ in range(4)]
        errors = [np.nan for _ in range(4)]

    if ERRORE == 0:
        y_fit_gaussian = negative_gaussian(wl, *popt_gaussian)
        errors = np.sqrt(np.diag(pcov))
        continuum_g=popt_gaussian[0]
        amplitude_g=popt_gaussian[1]
        line_center_g=popt_gaussian[2]
        sigma=popt_gaussian[3]
        continuum_g_error = errors[0]
        amplitude_g_error = errors[1]
        line_center_g_error = errors[2]
        sigma_error = errors[3]

    # Fit Lorentzian
    ERRORE=0
    try:
        popt_lorentzian, pcov = curve_fit(negative_lorentzian, wl, flux, p0=[const_guess, amplitude_guess, mean_guess, stddev_guess], maxfev=2000)
    except RuntimeError as e:
        ERORRE = 1
        popt_lorentzian= [np.nan for _ in range(4)]
        pcov = [[np.nan for _ in range(4)] for _ in range(4)]
        errors = [np.nan for _ in range(4)]
        delta_wl=np.nan
        rv=np.nan

    if ERRORE == 0:
        y_fit_lorentzian = negative_lorentzian(wl, *popt_lorentzian)
        errors = np.sqrt(np.diag(pcov))
        continuum=popt_lorentzian[0]
        amplitude=popt_lorentzian[1]
        line_center=popt_lorentzian[2]
        gamma=popt_lorentzian[3]
        continuum_error = errors[0]
        amplitude_error = errors[1]
        line_center_error = errors[2]
        gamma_error = errors[3]
        delta_wl=line_center0-popt_lorentzian[2]
        rv = (line_center/line_center0 - 1)*const.c.to(u.km/u.s).value

        mask1 = (wl > (line_center - 2 * gamma)) & (wl < (line_center + 2 * gamma))
        wl_subset=wl[mask1]
        flux_subset = flux[mask1]/continuum
        y_fit_lorentzian_subset = y_fit_lorentzian[mask1]/continuum
        residuals_fit = flux_subset - y_fit_lorentzian_subset
        rms_residuals = np.std(residuals_fit)
        snr=abs(amplitude/continuum)/rms_residuals
        chisq = np.sum((residuals_fit / np.std(flux_subset)) ** 2)
        dof = len(wl) - len(popt_lorentzian)
        reduced_chisq = chisq / dof
        ss_res = np.sum(residuals_fit ** 2)
        ss_tot = np.sum((flux_subset - np.mean(flux_subset)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        flux1=1-flux_subset
#        flux1=flux_subset
        sum_flux = sum(flux1)
        npt_subset = len(flux_subset)
        skewness = sum(flux1[i] * (wl_subset[i] - line_center)**3 for i in range(npt_subset)) / (sum_flux * sigma**3)
        kurtosis = sum(flux1[i] * (wl_subset[i] - line_center)**4 for i in range(npt_subset)) / (sum_flux * sigma**4) - 3

        log=True
        if log == True:
            print('continuum         ',continuum, continuum_error)
            print('intensity         ',amplitude, amplitude_error)
            print('central wavelength',line_center,line_center_error)
            print('gamma   = width   ',gamma, gamma_error)
            print('sigma   = width   ',sigma, sigma_error)
            print('skewness          ',skewness)
            print('kurtosis          ',kurtosis)
            print('rms residuals     ',rms_residuals)
            print('SNR               ',snr)
            print('chi2              ',reduced_chisq)
            print('R2                ',r_squared)
            print('delta wavelength  ',delta_wl)
            print('RV                ',rv)
            
        plot = True
        if plot == True:
            # Use these indices to extract the corresponding wavelengths and fluxes
            peak_wavelengths = wl[peak_indices_i]
            peak_fluxes = flux[peak_indices_i] / continuum
            x1 = line_center0 - wl_span
            x2 = line_center0 + wl_span
            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
            manager = plt.get_current_fig_manager()
            manager.window.wm_geometry("1200x600+0+0")
            ax1.plot(wl,flux/continuum, label='Data')
            #ax1.plot(wl,1-smoothed_flux/continuum, label='Data')
            #ax1.plot(wl,1-smoothed_flux0/continuum, label='Data')
            ax1.plot(wl, 1-y_fit_lorentzian/continuum, label='Lorentzian Fit', linestyle='dashed')
            ax1.axvline(line_center0, color='black', linestyle='dashed', label='Center')
            ax1.axvline(line_center-gamma, color='red', linestyle='dashed', label='Center')
            ax1.axvline(line_center+gamma, color='red', linestyle='dashed', label='Center')
            ax1.plot(peak_wavelengths, 1-peak_fluxes, 'o', label='Peaks')
            ax1.set_xlim(x1, x2)
            ax1.set_ylabel('Flux')
            ax1.grid(True)
            ax2.plot(wl, residuals_fit/continuum, label='Residuals')
            ax2.plot(wl, (flux-y_fit_lorentzian)/continuum, label='Residuals')
            ax2.axvline(line_center0, color='black', linestyle='dashed', label='Center')
            ax2.axvline(line_center-gamma, color='red', linestyle='dashed', label='Center')
            ax2.axvline(line_center+gamma, color='red', linestyle='dashed', label='Center')
            ax2.set_xlim(x1, x2)
            ax2.set_xlabel('Wavelength [AA]')
            ax2.set_ylabel('Residuals')
            ax2.grid(True)
            plt.tight_layout()
            plt.show()
    else:
        figatomare = True

    return rv,reduced_chisq,gamma,skewness,kurtosis
# ------------------------------  END   FUNCTION ------------------------------ #



# ------------------------------  BEGIN FUNCTION ------------------------------ #
def convert_to_air_wavelengths(wavelengths_vac):
    #Convert vacuum wavelengths to air wavelengths using the formula
    #from Donald Morton (2000, ApJ. Suppl., 130, 403)
    #n = 1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2), where s = 104 / λvac and λvac is in Ångströms.
    #The conversion is then: λair = λvac / n.
    #see https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
    
    # Calculate s squared
    s = 10**4 / wavelengths_vac
    s2 = s**2
    # Calculate the refractive index 
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2)
    # Convert to air wavelengths
    wavelengths_air = wavelengths_vac / n
    
    return wavelengths_air
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def negative_gaussian(x, const, amplitude, mean, stddev):
    gauss = const - amplitude * np.exp(-((x - mean) / stddev) ** 2)
    return gauss
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def negative_lorentzian(x, const, amplitude, x0, gamma):
    lorentz = const-amplitude * (gamma / ((x - x0)**2 + gamma**2))
    return lorentz
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def fit_absorption_line(wavelengths, intensities, wc,ws, flagd):
    plot=flagd
    log=flagd
    line_center0 = wc
    wl_span = ws
    mask = (wavelengths > (line_center0 - wl_span)) & (wavelengths < (line_center0 + wl_span))
    wl=wavelengths[mask]
    flux = intensities[mask]
    npt=len(wavelengths[mask])
    
    #for attempt in range(5):
    const_guess = np.mean(flux)
    amplitude_guess = np.min(flux) - const_guess
    mean_guess = line_center0 #*(attempt+1)
    stddev_guess = 0.1

    #fit Gaussian
    ERRORE=0
    try:
        popt_gaussian, pcov = curve_fit(negative_gaussian, wl, flux, p0=[const_guess, amplitude_guess, mean_guess, stddev_guess], maxfev=2000)
    except RuntimeError as e:
        ERORRE = 1
        popt_gaussian= [np.nan for _ in range(4)]
        pcov = [[np.nan for _ in range(4)] for _ in range(4)]
        errors = [np.nan for _ in range(4)]

    if ERRORE == 0:
        y_fit_gaussian = negative_gaussian(wl, *popt_gaussian)
        errors = np.sqrt(np.diag(pcov))
        sigma_error = errors[3]
        sigma=popt_gaussian[3]

    # Fit Lorentzian
    ERRORE=0
    try:
        popt_lorentzian, pcov = curve_fit(negative_lorentzian, wl, flux, p0=[const_guess, amplitude_guess, mean_guess, stddev_guess], maxfev=2000)
    except RuntimeError as e:
        ERORRE = 1
        popt_lorentzian= [np.nan for _ in range(4)]
        pcov = [[np.nan for _ in range(4)] for _ in range(4)]
        errors = [np.nan for _ in range(4)]
        delta_wl=np.nan
        rv=np.nan

    if ERRORE == 0:
        y_fit_lorentzian = negative_lorentzian(wl, *popt_lorentzian)
        errors = np.sqrt(np.diag(pcov))  # Extract the errors from the covariance matrix
    
        continuum=popt_lorentzian[0]
        amplitude=popt_lorentzian[1]
        line_center=popt_lorentzian[2]
        gamma=popt_lorentzian[3]
        continuum_error = errors[0]
        amplitude_error = errors[1]
        line_center_error = errors[2]
        gamma_error = errors[3]
        delta_wl=line_center0-popt_lorentzian[2]
        rv = (line_center/line_center0 - 1)*const.c.to(u.km/u.s).value
        mask1 = (wl > (line_center - 2 * gamma)) & (wl < (line_center + 2 * gamma))
        
        # Subset of flux and y_fit
        wl_subset=wl[mask1]
        flux_subset = flux[mask1]/continuum
        y_fit_lorentzian_subset = y_fit_lorentzian[mask1]/continuum
        
        residuals = flux_subset - y_fit_lorentzian_subset
        rms_residuals = np.std(residuals)
        snr=abs(amplitude/continuum)/rms_residuals
        
        chisq = np.sum((residuals / np.std(flux_subset)) ** 2)
        dof = len(wl) - len(popt_lorentzian)
        reduced_chisq = chisq / dof
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((flux_subset - np.mean(flux_subset)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        if r_squared<0 :
            print("************************************")
            print("************************************")
            print("* W A R N I N G       R2<0          ")
            print("************************************")
            print("************************************")
            #    break
            
    
    #if ERRORE == 0:
        #print('attempt ',attempt)
        #if attempt==1000:
        #        plot=False
        #        log=False
                
        flux1=1-flux_subset
        sum_flux = sum(flux1)
        npt_subset = len(flux_subset)
        
        skewness = sum(flux1[i] * (wl_subset[i] - line_center)**3 for i in range(npt_subset)) / (sum_flux * sigma**3)
        kurtosis = sum(flux1[i] * (wl_subset[i] - line_center)**4 for i in range(npt_subset)) / (sum_flux * sigma**4) - 3

        #if amplitude>0 and line_center0>line_center-gamma and line_center0<line_center+gamma and abs(line_center-line_center0)<0.25:
        if plot == True:    
    #       log=True
            if log == True:
                print('continuum         ',continuum, continuum_error)
                print('intensity         ',amplitude, amplitude_error)
                print('central wavelength',line_center,line_center_error)
                print('gamma   = width   ',gamma, gamma_error)
                print('sigma   = width   ',sigma)
                print('skewness          ',skewness)
                print('kurtosis          ',kurtosis)
                print('rms residuals     ',rms_residuals)
                print('SNR               ',snr)
                print('chi2              ',reduced_chisq)
                print('R2                ',r_squared)
                print('delta wavelength  ',delta_wl)
                print('RV                ',rv)
            #Blue: '#1F77B4'
            #Orange: '#FF7F0E'
            #Green: '#2CA02C'
            #Red: '#D62728'
            #Purple: '#9467BD'
    #        plot = True
     
            if plot == True:
                print(rv/const.c.to(u.km/u.s).value)
                x1 = line_center0 - wl_span
                x2 = line_center0 + wl_span
                #x1 = line_center0 - 2
                #x2 = line_center0 + 2
                lambda_values = absorption_lines_data['lambda_air']
                loggf_values = absorption_lines_data['log_gf']
                elements = absorption_lines_data['specie']
                mask_lines = (lambda_values > (line_center0 - wl_span)) & (lambda_values < (line_center0 + wl_span)) 
                #mask_lines2 = (lambda_values > (line_center0 - wl_span)) & (lambda_values < (line_center0 + wl_span)) & np.isnan(loggf_values)
                #new_wl = wl*(1+rv/const.c.to(u.km/u.s).value)
                #flux_interpolator = interp1d(new_wl, flux, kind='cubic', bounds_error=False, fill_value=0)
                #new_flux = flux_interpolator(wl)
                new_wl = wl*(1-rv/const.c.to(u.km/u.s).value)
                new_flux = flux
                print("Radial velocity (km/s):", rv)
                print("Speed of light (km/s):", const.c.to(u.km/u.s).value)
                print("Expected shift in Angstroms:", 6562 * (rv / const.c.to(u.km/u.s).value))
                print('mean shift in wavelength',np.mean(wl-new_wl),'mean difference in flux',np.mean(flux-new_flux))
                #for w, f, new_w, new_f in zip(wl, flux, new_wl, new_flux):
                #    print(w, f, new_w, new_f)
                print()
                num_plot = 1
                if num_plot == 1:
                    fig, ax1 = plt.subplots(1, 1)
                    manager = plt.get_current_fig_manager()
                    manager.window.wm_geometry("1400x1200+0+0")
                    original = False
                    if original == True:
                        ax1.plot(new_wl,new_flux/continuum, color='orange', linestyle ='dashed', label='Data')
                        ax1.axvline(line_center0, color='#00FF00', linestyle='dashed', label='Center')
                        ax1.axvline(line_center , color='#FF00FF', linestyle='dashed', label='Center')
                        ax1.axvline(line_center-gamma, color='#FF00FF', linestyle='dashed', label='Center')
                        ax1.axvline(line_center+gamma, color='#FF00CC', linestyle='dashed', label='Center')
                    else:
                        ax1.plot(new_wl,new_flux/continuum, label='Data')
                        ax1.plot(new_wl, y_fit_lorentzian/continuum, label='Lorentzian Fit', linestyle='dashed')
                        ax1.axvline(line_center0, color='#00FF00', linestyle='dashed', label='Center')
                        jj=-1
                        for lam, elem in zip(lambda_values[mask_lines], elements[mask_lines]):
                            jj= jj+1
                            ax1.axvline(lam, color = '#bbbbff', linestyle='dotted')
                            ax1.text(lam, ax1.get_ylim()[0]+0.0001*jj, elem, rotation=45, va='bottom')
                        #for lam, elem in zip(lambda_values[mask_lines2], elements[mask_lines2]):
                        #    ax1.axvline(lam, linestyle='dotted')
                        #    ax1.text(lam, ax1.get_ylim()[0], elem, rotation=45, va='bottom')
                    ax1.set_xlim(x1, x2)
                    ax1.set_ylabel('Flux')
                    ax1.grid(True)
                    plt.tight_layout()
                    plt.show()
                elif num_plot == 2:
                    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
                    manager = plt.get_current_fig_manager()
                    manager.window.wm_geometry("1800x600+0+0")
                    original = False
                    if original == True:
                        ax1.plot(new_wl,new_flux/continuum, color='orange', linestyle ='dashed', label='Data')
                        ax1.axvline(line_center0, color='#00FF00', linestyle='dashed', label='Center')
                        ax1.axvline(line_center , color='#FF00FF', linestyle='dashed', label='Center')
                        ax1.axvline(line_center-gamma, color='#FF00FF', linestyle='dashed', label='Center')
                        ax1.axvline(line_center+gamma, color='#FF00CC', linestyle='dashed', label='Center')
                    else:
    #                    ax1.plot(wl,flux/continuum, color ='#fdfdfd', label='Data')
                        ax1.plot(new_wl,new_flux/continuum, label='Data')
                        ax1.plot(new_wl, y_fit_lorentzian/continuum, label='Lorentzian Fit', linestyle='dashed')
                        ax1.axvline(line_center0, color='#00FF00', linestyle='dashed', label='Center')
                        for lam, elem in zip(lambda_values[mask_lines], elements[mask_lines]):
                            ax1.axvline(lam, color='#00e000', linestyle='dotted')
                            ax1.text(lam, ax1.get_ylim()[0], elem, rotation=45, va='top', color='#00a000')
                        for lam, elem in zip(lambda_values[mask_lines2], elements[mask_lines2]):
                            ax1.axvline(lam, color='#d000b0', linestyle='dotted')
                            ax1.text(lam, ax1.get_ylim()[0], elem, rotation=45, va='bottom', color='#8000b0')
                    ax1.set_xlim(x1, x2)
                    ax1.set_ylabel('Flux')
                    #ax1.set_title('H-alpha Line')
                    ax1.grid(True)
                    #ax1.legend()
                    ax2.plot(new_wl, (flux - y_fit_lorentzian)/continuum, label='Residuals')
                    ax2.axvline(line_center0, color='#00FF00', linestyle='dashed', label='Center')
                    #ax2.axvline(line_center-gamma, color='red', linestyle='dashed', label='Center')
                    #ax2.axvline(line_center+gamma, color='red', linestyle='dashed', label='Center')
                    for lam, elem in zip(lambda_values[mask_lines], elements[mask_lines]):
                        ax2.axvline(lam, color='#000080', linestyle='dotted')
                    ax2.set_xlim(x1, x2)
                    ax2.set_xlabel('Wavelength [AA]')
                    ax2.set_ylabel('Residuals')
                    ax2.grid(True)
                    plt.tight_layout()
                    plt.show()
                else:
                    figatomare = True
            
        new_wl = wl*(1-rv/const.c.to(u.km/u.s).value)
        new_flux = flux
        speck=new_flux/continuum
        speck_wl=new_wl
    else:
        figatomare = True

    return rv,reduced_chisq,gamma,skewness,kurtosis,speck,speck_wl
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def true_azimuth(angle):
    # convert ESO azimuth convention to astronomical definition of azimuth
    # till 2006.75 ESO used azimuth ranging from 180 to 540 clockwise (N=180, E=270, S=360, W=450)
    # from 2006.75 ESO used azimuth ranging from   0 to 360 clockwise (N=180, E=270, S=  0, W= 90)
    # the correction below (N=0, E=90, S=180, W=270)
    true_azimuth_angle = (angle + 180) % 360
    return true_azimuth_angle
# ------------------------------  END   FUNCTION ------------------------------ #

 
# ------------------------------  BEGIN FUNCTION ------------------------------ #
def convert_from_pseudosexagesimal_to_degrees(hhmmss, ddmmss):
# Example
#hhmmss = 123456.43234
#ddmmss = -1234.5678
#(hhmmss_str, hhmmss_deg), (ddmmss_str, ddmmss_deg) = convert_to_sexagesimal_and_degrees(hhmmss, ddmmss)
#print(f"HHMMSS Sexagesimal Representation: {hhmmss_str}")  
#print(f"HHMMSS Degrees: {hhmmss_deg}")
#print(f"DDMMSS Sexagesimal Representation: {ddmmss_str}")  
#print(f"DDMMSS Degrees: {ddmmss_deg}")

    # Function to convert a single part
    def convert_part(number, is_degrees=False):
        # Separating the integer and fractional parts
        integer_part = int(abs(number))
        fractional_part = abs(number) - integer_part

        # Converting integer part
        hours_or_degrees = integer_part // 10000
        minutes = (integer_part // 100) % 100
        seconds = integer_part % 100 + fractional_part
        
        # Forming the sexagesimal string representation
        sexagesimal_str = f"{int(hours_or_degrees):02d}:{int(minutes):02d}:{seconds:02f}".rstrip('0').rstrip('.')
        if number < 0 and is_degrees:
            sexagesimal_str = '-' + sexagesimal_str
        
        # Convert to degrees
        total_seconds = hours_or_degrees * 3600 + minutes * 60 + seconds
        if is_degrees:
            degrees = total_seconds / 3600
            if number < 0:
                degrees = -degrees
        else:
            degrees = total_seconds * 15 / 3600
        
        return sexagesimal_str, degrees
    
    # Convert hh:mm:ss.xs
    sexagesimal_hhmmss, degrees_hhmmss = convert_part(hhmmss)
    
    # Convert dd:mm:ss.xs
    sexagesimal_ddmmss, degrees_ddmmss = convert_part(ddmmss, True)
    
    return degrees_hhmmss, degrees_ddmmss,  sexagesimal_hhmmss, sexagesimal_ddmmss
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def testfilefits(myfile):
    if os.path.exists(myfile) and os.path.getsize(myfile) > 0:
        with open(myfile, 'rb') as f:
            first_line = f.read(6) # read the first 6 characters
            if first_line.startswith(b'SIMPLE'):
                isfits = True
            else:
                print(f'{myfile} is_not_a_valid_FITS_file')
                isfits = False
    else:
        print(f'ERR {myfile.replace("/scratch/mbarbier/mylocal/machines/Desktop_media/SeagateHub/DATA/HARPS/harpstar/data/","")} not_found_or_empty')
        isfits = False
    return isfits
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def calculate_star_altaz(ra, dec, geolon, geolat, elevation, mjd, humidity=None, temperature=None, pressure=None):
    star_coords = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    obs_location = EarthLocation(lon=geolon * u.deg, lat=geolat * u.deg, height=elevation * u.m)
    obs_time = Time(mjd, format='mjd')
    altaz_frame = AltAz(obstime=obs_time, location=obs_location)

    star_altaz = star_coords.transform_to(altaz_frame)
    star_azimuth0 = star_altaz.az
    star_altitude0 = star_altaz.alt

    altaz_frame_with_atmosphere = AltAz(obstime=obs_time, location=obs_location, obswl=550 * u.nm,temperature=temperature * u.deg_C, pressure=pressure * u.mbar, relative_humidity=humidity)
    star_altaz = star_coords.transform_to(altaz_frame_with_atmosphere)
    star_azimuth = star_altaz.az
    star_altitude = star_altaz.alt

    return star_azimuth, star_altitude, star_azimuth0, star_altitude0
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def calculate_sun_position(mjd, lon, lat, elevation):
    location = EarthLocation(lon=lon*u.deg, lat=lat*u.deg, height=elevation*u.m)
    time = Time(mjd, format='mjd')

    sun_altaz = get_sun(time).transform_to(AltAz(obstime=time, location=location))
    sun_ra_dec = get_sun(time).transform_to('icrs')

    sun_altitude = sun_altaz.alt
    sun_azimuth = sun_altaz.az

    return sun_ra_dec.ra, sun_ra_dec.dec, sun_altitude, sun_azimuth
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def calculate_moon_position(mjd, lon, lat, elevation):
    location = EarthLocation(lon=lon*u.deg, lat=lat*u.deg, height=elevation*u.m)
    time = Time(mjd, format='mjd')

    moon_altaz = get_moon(time).transform_to(AltAz(obstime=time, location=location))
    moon_ra_dec = get_moon(time).transform_to('icrs')

    moon_altitude = moon_altaz.alt
    moon_azimuth = moon_altaz.az

    return moon_ra_dec.ra, moon_ra_dec.dec, moon_altitude, moon_azimuth
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def butter_lowpass(cutoff, nyq_freq, order=4):
    normal_cutoff = float(cutoff) / nyq_freq
    b, a = signal.butter(order, normal_cutoff, btype='lowpass')
    return b, a
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def butter_lowpass_filter(data, cutoff_freq, nyq_freq, order=4):
    b, a = butter_lowpass(cutoff_freq, nyq_freq, order=order)
    y = signal.filtfilt(b, a, data)
    return y
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def butter_highpass(cutoff, nyq_freq, order=4):
    normal_cutoff = float(cutoff) / nyq_freq
    b, a = signal.butter(order, normal_cutoff, btype='highpass')
    return b, a
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def butter_highpass_filter(data, cutoff_freq, nyq_freq, order=4):
    b, a = butter_highpass(cutoff_freq, nyq_freq, order=order)
    y = signal.filtfilt(b, a, data)
    return y
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def statistics(y):
    max_y = np.nanmax(y)
    min_y = np.nanmin(y)
    sum_y = np.nansum(y)
    mean_y = np.nanmean(y)
    median_y = np.nanmedian(y)
    rms_y = np.nanstd(y, ddof=1)
    skew_y = stats.skew(y / max(y), nan_policy='omit')
    kurt_y = stats.kurtosis(y / max(y), nan_policy='omit')
    return max_y, min_y, sum_y, mean_y, median_y, rms_y, skew_y, kurt_y
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def process_spectrum(flux, wl, j1, j2):
    x = flux[0][j1:j2]
    wl_part = wl[0][j1:j2]
    max_value, min_value, sum_value, mean_value, median_value, rms_value, skew_value, kurt_value = statistics(x)
    coarse_resolution = 100.0  # Angstrom
    nnn = int(round((wl_part[-1] - wl_part[0]) / coarse_resolution))
    nwl = np.linspace(wl_part[0], wl_part[-1], nnn)
    spectrum_LR = np.zeros_like(nwl)
    for i in range(nnn):
        start_idx = np.searchsorted(wl_part, nwl[i] - coarse_resolution / 2)
        end_idx = np.searchsorted(wl_part, nwl[i] + coarse_resolution / 2)
        spectrum_LR[i] = np.sum(x[start_idx:end_idx]) / (end_idx - start_idx)
    f = interpolate.interp1d(nwl, spectrum_LR, kind='linear')
    spectrum_envelope = f(wl_part)
    spectrum_flat_0 = x / spectrum_envelope
    spec_median = np.nanmedian(spectrum_flat_0)
    spectrum_flat = spectrum_flat_0 / spec_median
    rms_f, min_f, sum_f, mean_f, median_f, rms_f, skew_f, kurt_f = statistics(spectrum_flat)
    return max_value, min_value, sum_value, mean_value, median_value, rms_value, skew_value, kurt_value, rms_f, min_f, sum_f, mean_f, median_f, rms_f, skew_f, kurt_f
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def process_file(file_list, file_type, datakw):
    spectral_masks = ["G2", "K0", "K5", "M2", "M4"]
    file_data = {}
    for keyword in datakw[file_type]['keywords']:
        column_name = keyword.replace('HIERARCH ESO ', '').replace(' ', '_').lower()
        for mask in spectral_masks:
            file_data[f'{column_name}_{mask}_{file_type}'] = None
    filefits = next((f for f in file_list if f"_{file_type}_" in f), None)
    if filefits:
        isfits = testfilefits(filefits)
        if isfits == True:
            hdu=fits.open(filefits, ignore_missing_end=True)
            #print(f'{file_type} {filefits.replace("/scratch/mbarbier/mylocal/machines/Desktop_media/SeagateHub/DATA/HARPS/harpstar/data/","")}')
            hdr0 = hdu[0].header
            for keyword in datakw[file_type]['keywords']:
                column_name = keyword.replace('HIERARCH ESO ', '').replace(' ', '_').lower()
                for mask in spectral_masks:
                    if mask in filefits:
                        try:
                            keyword_value = hdr0.get(keyword, None)
                            keyword_comment = hdr0.comments[keyword]
                            file_data[f'{column_name}_{mask}_{file_type}'] = {'value': keyword_value,'comment': keyword_comment}
                        except KeyError:
                            continue
            hdu.close()
    return file_data
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def process_gui_file(gui_files, datakw):
    gui_data = {}

    # Initialize gui_data dictionary
    for keyword in datakw['gui']['keywords']:
        column_name = keyword.replace('HIERARCH ESO ', '').replace(' ', '_').lower()
        gui_data[column_name] = {'value': None, 'comment': None}

    # If gui_files is not empty
    if gui_files:
        gui_file = gui_files[0]
        isfits = testfilefits(gui_file)
        if isfits:
            # Open the file
            with fits.open(gui_file) as harps_list_obs:
                hdr0 = harps_list_obs[0].header
                # Loop through the keywords
                for keyword in datakw['gui']['keywords']:
                    column_name = keyword.replace('HIERARCH ESO ', '').replace(' ', '_').lower()
                    try:
                        # Try to access the keyword
                        keyword_value = hdr0.get(keyword, None)
                        keyword_comment = hdr0.comments[keyword]
                        gui_data[column_name] = {'value': keyword_value, 'comment': keyword_comment}
                    except KeyError:
                        # Keyword not found, skip to the next one
                        print(f"Keyword '{keyword}' not found, skipping.")
                        continue

    return gui_data
# ------------------------------  END   FUNCTION ------------------------------ #



#********************************************************************************
#
#  CONFIGURATION
#
#********************************************************************************


# Location to store the downloaded files
storage_path = '/scratch/mbarbier/mylocal/machines/Desktop_media/SeagateHub/DATA/HARPS/harpstar/data'
cache_path = os.path.join(storage_path, '.cache')
exec_path = '/home/mbarbier/mylocal/harpsrvcatalog/output'


# Minimum and maximum usable wavelengths in the blue and red CCDs in Angstrom
wl_min_b = 3783.0
wl_max_b = 5304.0
bw_b = wl_max_b - wl_min_b  # 1521.0
wl_min_r = 5338.0
wl_max_r = 6912.0
bw_r = wl_max_r - wl_min_r  # 1574.0

file_names = ['adp', 'ccf', 'bis', 'gui']
base_filename = '{}_small.hdr'
spectral_masks = ["G2", "K0", "K5", "M2", "M4"]


#
stellar_templates = pd.read_csv("list.phoenix_stellar_templates.csv")
w1 = 3781.05
w2 = 6913.87
tpl_wl = np.arange(w1, w2, 0.01)
for index, row in stellar_templates.iterrows():
    spt = row['spt']
    template_file = row['template_file']
    
    with fits.open(template_file) as harps_list_obs:
        spectpl = harps_list_obs[1].data
        tpl_wl_vacuum = spectpl.field('WAVE')
        tpl_wl_air = convert_to_air_wavelengths(tpl_wl_vacuum)
        flux = spectpl.field('FLUX')
        interpolated_flux = np.interp(tpl_wl, tpl_wl_air, flux)
        if 'M0' in spt:
            tpl_M = interpolated_flux
        if 'K2' in spt:
            tpl_K = interpolated_flux
        if 'G2' in spt:
            tpl_G = interpolated_flux
        if 'F5' in spt:
            tpl_F = interpolated_flux
        if 'A5' in spt:
            tpl_A = interpolated_flux
        if 'B8' in spt:
            tpl_B = interpolated_flux
templates_spectra = [tpl_M, tpl_K, tpl_G, tpl_F, tpl_A, tpl_B]



#

#with fits.open('spectral_line_list_coluzzi_harps.fits') as absorption_lines:
#    absorption_lines_data = absorption_lines[1].data

with fits.open('nist_lines_output.fits') as absorption_lines:
    absorption_lines_data = absorption_lines[1].data
    


# Dictionary to store keywords and comments for different file types
datakw = {}

# Loop through each file type
for file_type in file_names:
    # Lists to store the keyword names and comments for the current file type
    keywords_list = []
    comments_list = []
    # Read the file line by line
    with open(base_filename.format(file_type), 'r') as f:
        for line in f:
            # Split the line based on '=' character
            keyword_part, value_comment_part = line.split('=', 1)
            # Strip whitespaces from keyword name and add to the list
            keywords_list.append(keyword_part.strip())
            # Split the second part based on '/' character to separate the value and comment
            value_part, comment_part = value_comment_part.split('/', 1)
            # Strip whitespaces from comment and add to the list
            comments_list.append(comment_part.strip())
    # Store the lists in the dictionary for the current file type
    datakw[file_type] = {'keywords': keywords_list, 'comments': comments_list}

# Open the FITS file
with fits.open('harps_simbad_crossmatch_nearest.fits') as harps_list_obs:
    harps_list_obs_data = harps_list_obs[1].data
    indices = np.arange(len(harps_list_obs_data))
    # Shuffle the indices
    np.random.shuffle(indices)
    # Use the shuffled indices to access the rows in random order
    shuffled_harps_list_obs_data = harps_list_obs_data[indices]

#print('begin masking')

#mask_spt = [not row['spt_obj'].startswith(('O','B','A','F', 'G', 'K', 'M', 'N', '-')) for row in shuffled_harps_list_obs_data]
#mask_spt = [row['spt_obj'].startswith(('C','W','Q','R', 'S', '?')) for row in shuffled_harps_list_obs_data]
#mask_acen = [(row['ra']>219.84) & (row['ra']<219.9) & (row['dec']<-60.826) & (row['dec']>-60.844)  for row in harps_list_obs_data]
#mask_bethyi= [(row['ra']>6.34) & (row['ra']<6.48) & (row['dec']>-77.265) & (row['dec']<-77.245)  for row in harps_list_obs_data]
#mask_pole= [(row['dec']>-70) & (row['dec']<-90)  for row in harps_list_obs_data]
#mask=mask_acen
#mask=mask_bethyi
#mask=mask_pole

## Apply the mask to select only the rows that meet the conditions
#selected_rows = [row for row, m in zip(harps_list_obs_data, mask) if m]
#selected_rows = [row for row, m in zip(shuffled_harps_list_obs_data, mask) if m]
## Sort the selected_rows by 'mjd_obs' in descending order
#selected_rows = sorted(selected_rows, key=lambda x: x['mjd_obs'], reverse=False)

#print(len(selected_rows))
#for row in selected_rows:
#    print(row['spt_obj'])
    
#print('end masking')
#sys.exit()
k = 0
speck = None
#speck_wl = None
# Loop over all rows

#for row in harps_list_obs_data:
#
#for row in selected_rows:
#
for row in shuffled_harps_list_obs_data:
    k += 1
    #print()
    #if(k>1000):
    #    sys.exit()
    dp_id_raw = row['dp_id_raw']
    archive_id_spectra0 = row['archive_id_spectra'].strip()
    archive_id_ancillary0 = row['archive_id_ancillary'].strip()

    if len(archive_id_spectra0) > 0:
        archive_id_spectra = archive_id_spectra0 + '.fits'
        dateobs = dp_id_raw[6:16]
        spectra_file = os.path.join(storage_path, dateobs, archive_id_spectra)
        all_results = {}
        #adp_data = {}
        #for keyword in datakw['adp']['keywords']:
        #    column_name = keyword.replace('HIERARCH ESO ', '').replace(' ', '_').lower()
        #    adp_data[column_name] = {'value': None,'comment': None}
        #adp_data['tel_targ_delta_epoch'] = {'value': None,'comment': None}
        #adp_data['tel_targ_alpha_epoch'] = {'value': None,'comment': None}
        #adp_data['mjd-obs_center'] = {'value': None,'comment': None}
        #adp_data['diff_ra_ep_pnt'] = {'value': None,'comment': None}
        #adp_data['diff_de_ep_pnt'] = {'value': None,'comment': None}
        #adp_data['diff_ra_ob_pnt'] = {'value': None,'comment': None}
        #adp_data['diff_de_ob_pnt'] = {'value': None,'comment': None}
        #adp_data['num_pmra'] = {'value': None,'comment': None}
        #adp_data['num_pmde'] = {'value': None,'comment': None}
        isfits = testfilefits(spectra_file)
        if isfits == True:
            harps_list_obs = fits.open(spectra_file)
            #print(f'adp {spectra_file.replace("/scratch/mbarbier/mylocal/machines/Desktop_media/SeagateHub/DATA/HARPS/harpstar/data/","")}')
            hdr0 = harps_list_obs[0].header
            hdr1 = harps_list_obs[1].header
            spec = harps_list_obs[1].data
            adp_data = {}
            keywords = datakw['adp']['keywords']
            for keyword in keywords:
                keyword_value = hdr0.get(keyword, None)
                keyword_comment = hdr0.comments[keyword]
                column_name = keyword.replace('HIERARCH ESO ', '').replace(' ', '_').lower()
                adp_data[column_name] = {'value': keyword_value,'comment': keyword_comment}
                #print(adp_data[column_name])
            adp_data['tel_targ_delta_epoch'] = {'value': None,'comment': None}
            adp_data['tel_targ_alpha_epoch'] = {'value': None,'comment': None}
            adp_data['mjd-obs_center'] = {'value': None,'comment': None}
            adp_data['diff_ra_ep_pnt'] = {'value': None,'comment': None}
            adp_data['diff_de_ep_pnt'] = {'value': None,'comment': None}
            adp_data['diff_ra_ob_pnt'] = {'value': None,'comment': None}
            adp_data['diff_de_ob_pnt'] = {'value': None,'comment': None}
            adp_data['num_pmra'] = {'value': None,'comment': None}
            adp_data['num_pmde'] = {'value': None,'comment': None}
            
            #adp_data[''] = {'value': None,'comment': None}
            #print(adp_data['tel_targ_alpha']['value'],adp_data['tel_targ_delta']['value'])
            tel_targ_alpha_degrees, tel_targ_delta_degrees, tel_targ_alpha_sexag, tel_targ_delta_sexag = convert_from_pseudosexagesimal_to_degrees(adp_data['tel_targ_alpha']['value'],adp_data['tel_targ_delta']['value'])
            adp_data['tel_targ_alpha']['value'] = tel_targ_alpha_degrees
            adp_data['tel_targ_delta']['value'] = tel_targ_delta_degrees

            ra1 , de1, _ , _ = convert_from_pseudosexagesimal_to_degrees(adp_data['ins_adc1_ra']['value'],adp_data['ins_adc1_dec']['value'])
            adp_data['ins_adc1_ra']['value']  = ra1
            adp_data['ins_adc1_dec']['value'] = de1

            ra2 , de2, _ , _ = convert_from_pseudosexagesimal_to_degrees(adp_data['ins_adc2_ra']['value'],adp_data['ins_adc2_dec']['value'])
            adp_data['ins_adc2_ra']['value']  = ra2
            adp_data['ins_adc2_dec']['value'] = de2
            
            #j2000=51544.5
            pmra=adp_data['tel_targ_pma']['value']
            pmde=adp_data['tel_targ_pmd']['value']
            mjd=adp_data['mjd-obs']['value']
            tel_targ_alpha_epoch = tel_targ_alpha_degrees + pmra/3600*(mjd-51544.5)/365.25
            tel_targ_delta_epoch = tel_targ_delta_degrees + pmde/3600*(mjd-51544.5)/365.25
            adp_data['tel_targ_alpha_epoch']['value'] = tel_targ_alpha_epoch
            adp_data['tel_targ_delta_epoch']['value'] = tel_targ_delta_epoch
            pm_tot_ob = np.sqrt(pmra**2+pmde**2)
            adp_data['mjd-obs_center']['value']=mjd+adp_data['exptime']['value']*adp_data['ins_det1_tmmean']['value']
            diff_ra_ep_pnt = (adp_data['tel_targ_alpha_epoch']['value']-adp_data['ra']['value'])*3600*np.cos(adp_data['tel_targ_delta']['value']*np.pi/180)
            diff_de_ep_pnt = (adp_data['tel_targ_delta_epoch']['value']-adp_data['dec']['value'])*3600
            diff_ra_ob_pnt = (adp_data['tel_targ_alpha']['value']-adp_data['ra']['value'])*3600*np.cos(adp_data['tel_targ_delta']['value']*np.pi/180)
            diff_de_ob_pnt = (adp_data['tel_targ_delta']['value']-adp_data['dec']['value'])*3600
            #num_pmra = diff_ra_ep_pnt/(pmra*(mjd-51544.5)/365.25)
            #num_pmde = diff_de_ep_pnt/(pmde*(mjd-51544.5)/365.25)
            lon=adp_data['tel_geolon']['value']
            lat=adp_data['tel_geolat']['value']
            ele=adp_data['tel_geoelev']['value']
            humidity=adp_data['tel_ambi_rhum']['value']
            temperature=adp_data['tel_ambi_temp']['value']
            pressure=adp_data['tel_ambi_pres_start']['value']
            true_azimuth_angle = true_azimuth(adp_data['tel_az']['value'])
            adp_data['tel_az']['value'] = true_azimuth_angle

            sun_ra , sun_dec , sun_alt , sun_az  = calculate_sun_position(mjd, lon, lat, ele)
            moon_ra, moon_dec, moon_alt, moon_az = calculate_moon_position(mjd, lon, lat, ele)
            sun_coord = SkyCoord(ra=sun_ra, dec=sun_dec, frame='icrs')
            moon_coord = SkyCoord(ra=moon_ra, dec=moon_dec, frame='icrs')
            object_coord = SkyCoord(ra=tel_targ_alpha_epoch*u.deg, dec=tel_targ_delta_epoch*u.deg, frame='icrs')
            
            angular_distance_obj_sun = object_coord.separation(sun_coord)
            angular_distance_obj_moon = object_coord.separation(moon_coord)
            
            #print(angular_distance_obj_sun.deg,angular_distance_obj_moon.deg,sun_alt.deg)

            azimuth, altitude, azimuth0, altitude0 = calculate_star_altaz(tel_targ_alpha_epoch, tel_targ_delta_epoch, lon, lat, ele, mjd, humidity, temperature, pressure)
            diff_az  = (azimuth.deg-adp_data['tel_az']['value'])*3600*np.cos(altitude*np.pi/180)
            diff_alt = (altitude.deg-adp_data['tel_alt']['value'])*3600


            all_results.update(adp_data)
            wl = spec.field('WAVE')
            flux = spec.field('FLUX')
            wl_ = spec.field('WAVE').flatten()
            flux_ = spec.field('FLUX').flatten()
            npoints = flux.size
            
            #print(adp_data['object']['value'])
            # Tecnezio
            #popt,deltawl,rvha = fit_absorption_line(wl_, flux_,4238.19,3.0)
            #popt,deltawl,rvha = fit_absorption_line(wl_, flux_,4262.27,3.0)
            #popt,deltawl,rvha = fit_absorption_line(wl_, flux_,4297.00,3.0)
            
            print()
            print(row['object'],row['spt_obj'])
            #rv,reduced_chisq,gamma,skewness,kurtosis = fit_absorption_line2(wl_, flux_,6562.817,30.0)
            
            rv,reduced_chisq,gamma,skewness,kurtosis,speck1,speck_wl = fit_absorption_line(wl_, flux_,6562.817,20.0, False)
            
            if k==1:
                speck = np.zeros_like(speck1)

            print(k)
            speck=((k-1)*speck+speck1)/k
            
            if np.mod(k,100)==0:
                fig, ax1 = plt.subplots(1, 1)
                ax1.plot(speck)
                plt.show()
            

            
            #new_wl = wl_*(1+rv/const.c.to(u.km/u.s).value)
            #flux_interpolator = interp1d(wl_, flux_, kind='cubic', bounds_error=False, fill_value=0)
            #new_flux = flux_interpolator(new_wl)
            #_,_,_,_,_ = fit_absorption_line(new_wl, new_flux,6562.817,10.0, True)

            
            #popt,deltawl,rvha = fit_absorption_line(wl_, flux_,4861.332,10.0)
            #continuum,intensity,centrawl,width1sig = popt
            #print(row['object'],",",row['mjd_obs'],",",rv,",",reduced_chisq,",",gamma,",",skewness,",",kurtosis)

            

            stat_spec_b = {
                'max_b': None, 'min_b': None, 'sum_b': None, 'mean_b': None, 'median_b': None, 'rms_b': None,
                'skew_b': None, 'kurt_b': None, 'max_bf': None, 'min_bf': None, 'sum_bf': None, 'mean_bf': None,
                'median_bf': None, 'rms_bf': None, 'skew_bf': None, 'kurt_bf': None
            }
            stat_spec_r = {
                'max_r': None, 'min_r': None, 'sum_r': None, 'mean_r': None, 'median_r': None, 'rms_r': None,
                'skew_r': None, 'kurt_r': None, 'max_rf': None, 'min_rf': None, 'sum_rf': None, 'mean_rf': None,
                'median_rf': None, 'rms_rf': None, 'skew_rf': None, 'kurt_rf': None
            }
            
            if npoints > 150000:
                j_values = [0, 0, 0, 0]
                ranges = [(0, 300), (152100, 152400), (155400, 155700), (npoints - 300, npoints - 1)]
                conditions = [(wl_min_b, 'b'), (wl_max_b, 'b'), (wl_min_r, 'r'), (wl_max_r, 'r')]
                for index, ((start, end), (wl_limit, band)) in enumerate(zip(ranges, conditions)):
                    for j in range(start, end):
                        wl_prev = wl[0][j - 1]
                        wl_curr = wl[0][j]
                        if wl_prev < wl_limit <= wl_curr:
                            j_values[index] = j
                            break
                j1, j2, j3, j4 = j_values
                
                if j1 > 0 and j2 > j1:
                    max_b, min_b, sum_b, mean_b, median_b, rms_b, skew_b, kurt_b, rms_bf, min_bf, sum_bf, mean_bf, median_bf, rms_bf, skew_bf, kurt_bf = process_spectrum(flux, wl, j1, j2)
                    stat_spec_b = {'max_b': max_b, 'min_b': min_b, 'sum_b': sum_b, 'mean_b': mean_b, 'median_b': median_b, 'rms_b': rms_b,'skew_b': skew_b, 'kurt_b': kurt_b, 
                                   'max_bf': rms_bf, 'min_bf': min_bf, 'sum_bf': sum_bf, 'mean_bf': mean_bf, 'median_bf': median_bf, 'rms_bf': rms_bf, 'skew_bf': skew_bf, 'kurt_bf': kurt_bf}

                if j3 > j2 and j3 < j4:
                    max_r, min_r, sum_r, mean_r, median_r, rms_r, skew_r, kurt_r, rms_rf, min_rf, sum_rf, mean_rf, median_rf, rms_rf, skew_rf, kurt_rf = process_spectrum(flux, wl, j3, j4)
                    stat_spec_r = {'max_r': max_r, 'min_r': min_r, 'sum_r': sum_r, 'mean_r': mean_r, 'median_r': median_r, 'rms_r': rms_r,'skew_r': skew_r, 'kurt_r': kurt_r, 
                                   'max_rf': rms_rf, 'min_rf': min_rf, 'sum_rf': sum_rf,'mean_rf': mean_rf, 'median_rf': median_rf, 'rms_rf': rms_rf, 'skew_rf': skew_rf, 'kurt_rf': kurt_rf}
            
            all_results.update(stat_spec_b)
            all_results.update(stat_spec_r)
            harps_list_obs.close()


######
######
######
######
######

        # Create patterns to search for "bis" and "ccf" files
        # Check if the length of archive_id_ancillary is greater than 0 and if there is at least one "bis" or one "ccf" file
        # Process CCF file of this type
        ccf_pattern = os.path.join(storage_path, dateobs, dp_id_raw + "_ccf*_A.fits")
        ccf_files = glob.glob(ccf_pattern)
        ccf_data = process_file(ccf_files, 'ccf', datakw)
        all_results.update(ccf_data)

        # Process BIS file of this type
        bis_pattern = os.path.join(storage_path, dateobs, dp_id_raw + "_bis*_A.fits")
        bis_files = glob.glob(bis_pattern)
        bis_data = process_file(bis_files, 'bis', datakw)
        all_results.update(bis_data)

        # Process GUIDE file
        gui_pattern = os.path.join(storage_path, dateobs, dp_id_raw + "_INT_GUIDE.fits")
        gui_files = glob.glob(gui_pattern)
        gui_data = process_gui_file(gui_files, datakw)
        all_results.update(gui_data)
        
        
        #print(all_results)
        #print(k)










sys.exit()
#---------------   END PROGRAM













#radvel = False
#if radvel == True:
    #em1=np.min(wl) 
    #em2=np.max(wl)
    #ix1=np.argmin(np.abs(tpl_wl-em1))
    #ix2=np.argmin(np.abs(tpl_wl-em2))
    #tpl_wl_cut=tpl_wl[ix1:ix2+1]

    #results = []
    #max_velocity = +1000  # km/s
    #min_velocity = -1000  # km/s
    #deltalamba_over_lambda = 0.01 / 5500.0
    #vel_c = const.c.to(u.km/u.s).value
    #kost = (deltalamba_over_lambda * vel_c)
    ##print(80/kost)
    ##sys.exit()
    #kk=0
    #for flux_template in templates_spectra:
        #tpl_flux_cut = flux_template[ix1:ix2+1]
        #kk += 1
        #if kk >1:
            #break
        ## Calculate the maximum allowed lag for the desired velocity range
        #max_lag = int(max_velocity / kost )
        #min_lag = int(min_velocity / kost )
        ##print(min_lag,max_lag)
        
        ## Limit the range of lags
        #lag = np.arange(min_lag, max_lag + 1)
        #velocities = lag * kost
        #lag = np.arange(-npoints + 1, npoints)

        #correlation = correlate(flux_ - np.mean(flux_), tpl_flux_cut - np.mean(tpl_flux_cut), mode='full')
        ## Since the correlation is calculated for a larger range of lags than needed,
        ## slice the correlation array to match the range of lags of interest.
        #start_index = npoints + min_lag
        #end_index = npoints + max_lag + 1
        #correlation = correlation[start_index:end_index]
        
        #radial_velocity = velocities[np.argmax(correlation)]
        #peaks, _ = find_peaks(correlation)
        #results_half = peak_widths(correlation, peaks, rel_height=0.5)
        #fwhm = results_half[0]
        ##print(fwhm)
        #results.extend([radial_velocity, fwhm[0]])
        #rvu=adp_data['tel_targ_radvel']['value']
        #rv_obs=radial_velocity-adp_data['drs_berv']['value']
        #berv=adp_data['drs_berv']['value']
        #rv=rv_obs+berv
        #print(adp_data['date-obs']['value'])
        #print('RVuser        = ',rvu)
        #print('RV            = ',rv)
        #print('RVu-BERV-RV   = ',rvu-rv)
        
        #figures = True
        #if figures == True:
            #plt.figure(figsize=(12, 6))
            #plt.subplot(2, 1, 1)
            #plt.plot(wl_, flux_/np.mean(flux_), label='flux')
            #plt.plot(tpl_wl_cut, tpl_flux_cut/np.mean(tpl_flux_cut)+2, label='tpl_flux_cut')
            #plt.xlabel('Wavelength')
            #plt.ylabel('Flux')
            #plt.legend()
            
            #plt.subplot(2, 1, 2)
            #plt.plot(velocities, correlation, label='correlation')
            #plt.xlabel('Velocity (km/s)')
            #plt.ylabel('Correlation')
            #plt.legend()
            
            #plt.show()
