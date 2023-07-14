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
from collections import defaultdict
import csv
import datetime
#
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#
from scipy import signal, interpolate, stats
import scipy.stats as st
from scipy.signal      import correlate, find_peaks, peak_widths, savgol_filter, peak_prominences
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.ndimage     import gaussian_filter,gaussian_filter1d, uniform_filter1d
from scipy.optimize    import curve_fit
from scipy.integrate   import simps,quad
from scipy.special     import erf,j1
from scipy.stats       import linregress
#
import astropy.units as u
from astropy.io import fits
from astropy import constants as const
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_body, get_sun, get_moon
from astropy.time import Time


# Filter the warnings
logging.basicConfig(filename='error.log', level=logging.ERROR)
warnings.filterwarnings('ignore', category=UserWarning, append=True)
warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)

#********************************************************************************
#
#  CONSTANTS
#
#********************************************************************************

solar_mass=const.M_sun.to(u.g).value # g
solar_radius=const.R_sun.to(u.cm).value #cm
solar_luminosity=3.828e33 # erg/cm**2/s
plank_constant = const.h.value  # Planck's constant (J.s)
ligth_speed = const.c.to(u.km/u.s).value  # Speed of light (km/s)
boltzmann_constant = const.k_B.value  # Boltzmann's constant (J/K)
sb_constant=const.sigma_sb.value # Sfefan Boltzmann constant (W.K**-4.m**-2)
const_gravity = const.G.value # (m**3 kg**-1 s**-2)
wien_constant = const.b_wien.value # Wien displacement constant (m/K)
au = const.au.to(u.cm).value # astronomical unit (cm)
pc = 3085677581491367300.0 # parsec (cm)
ha_wl = 6562.817 # wavelegth of Halpha in AA

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
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def analyze_curvature(x, y):

    try:
        a, b, c = np.polyfit(x, y, 2)
    except np.linalg.LinAlgError:
        a = np.nan
        b = np.nan
        c = np.nan

    if a>0:
        curvature = 1
    elif a<0:
        curvature = -1
    else:
        curvature = 0
    return curvature
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def calc_atm_dispersion(altitude,temperature,pressure0,relhum0,year):
    pressure=pressure0*100 # from mbar to pascal
    relhum=relhum0/100 # from % to [0-1]

    wl_b = 0.3783  # µm
    wl_r = 0.6912  # µm
    
    co2=co2value(year)
    # Refractive indexes
    n_b = refractive_index(wl_b, temperature, pressure, relhum, co2)
    n_r = refractive_index(wl_r, temperature, pressure, relhum, co2)
    zenith_distance = np.radians(90 - altitude.deg)
    #delta is in arcsec
    delta=(n_b-n_r)*np.tan(zenith_distance)*180/np.pi*3600
    return delta
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def refractive_index(lambda_, temp_c, p, h, xc):
    # Ciddor formula
    # https://ui.adsabs.harvard.edu/abs/1996ApOpt..35.1566C/abstract
    # http://walter.bislins.ch/blog/media/Refractive%20index%20of%20air$2C%20new%20equations%20for%20the%20visible%20and%20near%20infrared.pdf
    # lambda_ = wavelength, 0.3 to 1.69 µm
    # t = temperature, -40 to +100 °C
    # p = pressure, 80000 to 120000 Pa
    # h = fractional humidity, 0 to 1
    # xc = CO2 concentration, 0 to 2000 ppm

    sigma = 1 / lambda_  # µm^-1
    sigma2 = sigma * sigma
    sigma4 = sigma2 * sigma2
    sigma6 = sigma4 * sigma2

    T = temp_c + 273.15  # Temperature °C -> K
    T2 = T**2
    t2_c = temp_c**2

    R = 8.314463  # gas constant, J/(mol·K)

    k0 = 238.0185  # µm^-2
    k1 = 5792105  # µm^-2
    k2 = 57.362  # µm^-2
    k3 = 167917  # µm^-2

    w0 = 295.235  # µm^-2
    w1 = 2.6422  # µm^-2
    w2 = -0.032380  # µm^-4
    w3 = 0.004028  # µm^-6

    A = 1.2378847e-5  # K^-2
    B = -1.9121316e-2  # K^-1
    C = 33.93711047
    D = -6.3431645e3  # K

    alpha = 1.00062
    beta = 3.14e-8  # Pa^-1,
    gamma = 5.6e-7  # °C^-2

    # saturation vapor pressure of water vapor in air at temperature T
    if temp_c >= 0:
        svp = math.exp(A * T2 + B * T + C + D / T)  # Pa
    else:
        svp = math.pow(10, -2663.5 / T + 12.537)

    # enhancement factor of water vapor in air
    f = alpha + beta * p + gamma * t2_c

    # molar fraction of water vapor in moist air
    xw = f * h * svp / p

    # refractive index of standard air at 15 °C, 101325 Pa, 0% humidity, 450 ppm CO2
    nas = 1 + (k1 / (k0 - sigma2) + k3 / (k2 - sigma2)) * 1e-8

    # refractive index of standard air at 15 °C, 101325 Pa, 0% humidity, xc ppm CO2
    naxs = 1 + (nas - 1) * (1 + 0.534e-6 * (xc - 450))

    # refractive index of water vapor at standard conditions (20 °C, 1333 Pa)
    nws = 1 + 1.022 * (w0 + w1 * sigma2 + w2 * sigma4 + w3 * sigma6) * 1e-8

    Ma = 1e-3 * (28.9635 + 12.011e-6 * (xc - 400))  # molar mass of dry air, kg/mol
    Mw = 0.018015  # molar mass of water vapor, kg/mol

    Za = compressibility(288.15, 101325, 0)  # compressibility of dry air
    Zw = compressibility(293.15, 1333, 1)  # compressibility of pure water vapor

    # Eq.4 with (T,P,xw) = (288.15, 101325, 0)
    rho_axs = 101325 * Ma / (Za * R * 288.15)  # density of standard air

    # Eq 4 with (T,P,xw) = (293.15, 1333, 1)
    rho_ws = 1333 * Mw / (Zw * R * 293.15)  # density of standard water vapor

    # two parts of Eq.4: rho = rho_a + rho_w
    Z = compressibility(T, p, xw)
    rho_a = p * Ma / (Z * R * T) * (1 - xw)  # rho_a: density of the dry component of the moist air
    rho_w = p * Mw / (Z * R * T) * xw  # rhow_w: density of the water vapor component

    nref = 1 + (rho_a / rho_axs) * (naxs - 1) + (rho_w / rho_ws) * (nws - 1)

    return nref
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def compressibility(T, p, xw):
    # T = temperature in K
    # p = pressure in Pa
    # xw = molar fraction of water vapor in moist air

    t = T - 273.15
    a0 = 1.58123e-6  # K·Pa^-1
    a1 = -2.9331e-8  # Pa^-1
    a2 = 1.1043e-10  # K^-1·Pa^-1
    b0 = 5.707e-6  # K·Pa^-1
    b1 = -2.051e-8  # Pa^-1
    c0 = 1.9898e-4  # K·Pa^-1
    c1 = -2.376e-6  # Pa^-1
    d = 1.83e-11  # K^2·Pa^-2
    e = -0.765e-8  # K^2·Pa^-2
    pt = p / T
    pt2 = pt * pt
    t2 = t * t
    xw2 = xw * xw
    compress = 1 - pt * (a0 + a1 * t + a2 * t2 + (b0 + b1 * t) * xw + (c0 + c1 * t) * xw2) + pt2 * (d + e * xw2)
    return compress
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def co2value(year):
    # return the CO2 value from the polynomial fit of the CO2 Keeling curve from Manua Kea
    # valid from 1958.2 to 2023.5
    # 
    # original data
    #https://raw.githubusercontent.com/sio-co2o2/keelingcurve_notebooks/main/preliminary_data/mlo/mlo_full_record_now_span.csv
    #
    # input:  decimal year
    # output: CO2 in ppm
    #
    # fit of original data
    # Given normalization constants
    #x_offset = 1958.2
    #x_scale = 65.27
    #y_offset = 313
    #y_scale = 112
    # fit coefficient
    #a1 = 0.537366984464985
    #b1 = 0.410621792713439
    #c1 = 0.0184694570646881
    # Calculate coefficients in original coordinates
    #a2 = (y_scale * a1) / (x_scale ** 2)
    #b2 = (y_scale * b1) / x_scale - (2 * y_scale * a1 * x_offset) / (x_scale ** 2)
    #c2 = (y_scale * a1 * x_offset ** 2) / (x_scale ** 2) - (y_scale * b1 * x_offset) / x_scale + y_scale * c1 + y_offset
    
    a2 = 0.014127384890899222
    b2 = -54.623884076959705
    c2 = 53107.43363687843
    
    co2 = a2*year**2+b2*year+c2
    if(year<1958.2 or year>2023.5):
        print("WARNING: CO2 values extrapolated outside valid range [1958.2-2023.5]")
    return co2
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def measure_equivalent_width(line_wavelength, wavelengths, flux, left_window, right_window):
    line_index = np.argmin(np.abs(wavelengths - line_wavelength))
    left_index = np.argmin(np.abs(wavelengths - (line_wavelength - left_window)))
    right_index = np.argmin(np.abs(wavelengths - (line_wavelength + right_window)))
    #continuum = np.mean(np.concatenate((flux[left_index:line_index], flux[line_index:right_index])))
    continuum = np.mean(flux[(line_index - left_index):(line_index + right_index)])
    delta_wavelength = np.mean(np.diff(wavelengths))
    ew = np.sum((1 - flux[left_index:right_index] / continuum) * delta_wavelength)
    return ew
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def linear_func(x, a, b):
    return a * x + b
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def quadratic_func(x, a, b, c):
    return a * x**2 + b *x + c
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #

def moving_avg(arr, window_size, polyorder):
    # Check window size is odd
    if window_size % 2 == 0:
        window_size += 1

    # Calculate moving averages using Savitzky-Golay filter
    moving_average = savgol_filter(arr, window_size, polyorder)  # 0th derivative for mean
    moving_squared_average = savgol_filter(arr**2, window_size, polyorder)
    moving_cubed_average = savgol_filter(arr**3, window_size, polyorder)
    moving_quad_average = savgol_filter(arr**4, window_size, polyorder)

    # Calculate moving variance (standard deviation)
    moving_variance = moving_squared_average - moving_average**2
    moving_std = np.sqrt(moving_variance)

    # Calculate moving skewness
    moving_skewness = (moving_cubed_average - 3*moving_average*moving_squared_average + 2*moving_average**3) / moving_std**3

    # Calculate moving kurtosis
    moving_kurtosis = ((moving_quad_average - 4*moving_average*moving_cubed_average + 6*moving_squared_average*moving_average**2 - 3*moving_average**4) / moving_std**4) - 3

    return moving_average, moving_std, moving_skewness, moving_kurtosis


# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def perform_fftold(wave, flux):
    N = len(flux)
    sampling_rate = wave[1] - wave[0]

    normalized_flux = flux - np.mean(flux)
    
    # FFT computation
    fft_values = np.fft.fft(normalized_flux)
    freq = np.fft.fftfreq(N, d=sampling_rate)
    
    # Selecting only positive frequencies
    positive_indices = np.where(freq > 0)
    fft_values = fft_values[positive_indices]
    freq = freq[positive_indices]
    
    # PSD computation
    psd = np.abs(fft_values) ** 2 / freq
    
    # Log-log scale
    log_psd = np.log10(psd)
    log_freq = np.log10(freq)
    
    # Smoothing PSD
    psd_smoothed = gaussian_filter1d(log_psd, sigma=30)
    
    # Fitting
    lf1 = -2
    lf2 = 1
    fit_indices = np.where((log_freq >= lf1) & (log_freq <= lf2))
    popt, pcov = curve_fit(linear_func, log_freq[fit_indices], psd_smoothed[fit_indices])
    fit_values = linear_func(log_freq[fit_indices], *popt)
    
    # Goodness of fit
    chi_squared = np.sum((psd_smoothed[fit_indices] - fit_values) ** 2 / fit_values)
    ss_res = np.sum((psd_smoothed[fit_indices] - fit_values) ** 2)
    ss_tot = np.sum((psd_smoothed[fit_indices] - np.mean(psd_smoothed[fit_indices])) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    perr = np.sqrt(np.diag(pcov))
    
    slope = popt[0]
    intercept = popt[1]
    slope_e = perr[0]
    intercept_e = perr[1]

    peaks, _ = find_peaks(np.abs(fft_values))
    # Get the index of the highest peak
    highest_peak_index = np.argmax(fft_values[peaks])
    # Retrieve the value and frequency of the highest peak
    highest_peak_value = fft_values[peaks][highest_peak_index]
    highest_peak_frequency = peaks[highest_peak_index]

    
    makefigure = False
    if makefigure == True:
        window_size = 15
        polord = 3
        moving_average, moving_std, moving_skew, moving_kurt = moving_avg(np.abs(fft_values), window_size ,polord)
    
        threshold = 3 
    
        # Find indices where the absolute FFT values are greater than the threshold times the standard deviation
        #peak_indices = np.where(np.abs(fft_values) > threshold * movig_std + moving_average)[0]
    
        #highest_peak_index = peaks[np.argmax(np.abs(fft_values[peaks]))]
        #highest_peak_frequency = freq[highest_peak_index]
        #highest_peak = np.abs(fft_values[highest_peak_index])
    
        #print(highest_peak_frequency,1/highest_peak_frequency)
        #print(peaks.size,fft_values.size)
        #print(np.mean(rolled_std)*100)
        #print(np.mean(rolled_mean))
        #print(np.mean(np.abs(fft_values)))
        #print(np.std(np.abs(fft_values)))
     
        # Plotting results
        plt.figure(figsize=(12, 8))

        plt.subplot(3,1,1)
        plt.plot(wave, normalized_flux)
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')

        plt.subplot(3,1,2)
        plt.xlim([0, 0.02])
        #plt.ylim([10, 20])
        plt.plot(freq, (np.abs(fft_values)))
        plt.plot(freq, moving_average)
        plt.plot(freq, moving_std)
        plt.plot(freq, moving_skew)
        plt.plot(freq, moving_kurt)
        plt.xlabel('Frequency')
        plt.ylabel('Power')

        plt.subplot(3,1,3)
        plt.loglog(freq, psd)
        plt.xlabel('Frequency')
        plt.ylabel('PSD')
 
        # Add the highest peak to the power and PSD plots
        #plt.subplot(3,1,2)
        #plt.plot(highest_peak_frequency, np.log10(highest_peak**2), 'ro')
        #plt.subplot(3,1,3)
        #plt.loglog(highest_peak_frequency, highest_peak**2, 'ro')

        # Add the corresponding sinusoid to the flux vs wavelength plot
        #plt.subplot(3,1,1)
        #sinusoid = np.sin(2 * np.pi * highest_peak_frequency * wave)
        #plt.plot(wave, sinusoid)
        plt.tight_layout()
        plt.show()
        plt.close()
    
    
    
    # Results
    return log_freq, log_psd, psd_smoothed, slope, intercept, slope_e, intercept_e, chi_squared, r_squared
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def perform_fft(wave, flux):
    try:
        N = len(flux)
#        if N <= 313000 or np.any(flux < 0):
        if N <= 310000:
            raise ValueError('Invalid data for FFT computation.')

        sampling_rate = wave[1] - wave[0]

        normalized_flux = flux - np.mean(flux)

        # FFT computation
        fft_values = np.fft.fft(normalized_flux)
        freq = np.fft.fftfreq(N, d=sampling_rate)

        # Selecting only positive frequencies
        positive_indices = np.where(freq > 0)
        fft_values = fft_values[positive_indices]
        freq = freq[positive_indices]

        # PSD computation
        psd = np.abs(fft_values) ** 2 / freq

        # Log-log scale
        log_psd = np.log10(psd)
        log_freq = np.log10(freq)

        # Smoothing PSD
        psd_smoothed = gaussian_filter1d(log_psd, sigma=30)

        # Fitting
        lf1 = -2
        lf2 = 1
        fit_indices = np.where((log_freq >= lf1) & (log_freq <= lf2))
        
        # Check for invalid values in the data for curve fitting
        if np.any(np.isnan(log_freq[fit_indices])) or np.any(np.isnan(psd_smoothed[fit_indices])):
            raise ValueError('Data for curve fitting contains NaNs or Infs.')

        popt, pcov = curve_fit(linear_func, log_freq[fit_indices], psd_smoothed[fit_indices])
        fit_values = linear_func(log_freq[fit_indices], *popt)

        # Goodness of fit
        chi_squared = np.sum((psd_smoothed[fit_indices] - fit_values) ** 2 / fit_values)
        ss_res = np.sum((psd_smoothed[fit_indices] - fit_values) ** 2)
        ss_tot = np.sum((psd_smoothed[fit_indices] - np.mean(psd_smoothed[fit_indices])) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        perr = np.sqrt(np.diag(pcov))

        slope = popt[0]
        intercept = popt[1]
        slope_e = perr[0]
        intercept_e = perr[1]

        peaks, _ = find_peaks(np.abs(fft_values))
        # Get the index of the highest peak
        highest_peak_index = np.argmax(fft_values[peaks])
        # Retrieve the value and frequency of the highest peak
        highest_peak_value = fft_values[peaks][highest_peak_index]
        highest_peak_frequency = peaks[highest_peak_index]
        highest_peak_wl = 1/highest_peak_frequency

        # Results
        return log_freq, log_psd, psd_smoothed, slope, intercept, slope_e, intercept_e, chi_squared, r_squared, highest_peak_value, highest_peak_wl

    except ValueError as e:
        return [np.nan]*11
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def clean_flux_and_calculate_magnitudes(flux, wavelengths, filters_data):
    # Clean flux: replace <0, NaN, Inf with 0
    cleaned_flux = np.where((flux < 0) | np.isnan(flux) | np.isinf(flux), 0, flux)
    
    # Calculate magnitude in each filter
    magnitudes = {}
    fluxes = {}
    for filter_name, filter_data in filters_data.items():
        filter_wavelength = filter_data['wavelength']
        filter_transmission = filter_data['transmission']
        
        transmission_interp = np.interp(wavelengths, filter_wavelength, filter_transmission)
        flux_filter = simps(cleaned_flux * transmission_interp, wavelengths)
        magnitude = -2.5 * np.log10(flux_filter)
        magnitudes[filter_name] = magnitude
        fluxes[filter_name] = flux_filter
        #print(filter_name,magnitude)
    
    return magnitudes,fluxes
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def calc_extcoeff(airmass):
    airmass_values = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9])
    airmass_columns = [f'airm{int(airmass * 10)}' for airmass in [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9]]
    transmissions = trasm_df[airmass_columns]
    spline_interpolator = interp1d(airmass_values, transmissions.values.T, kind='cubic', axis=0)
    spline_interpolated_transmissions = spline_interpolator(airmass)
    extcoeff=-np.log(spline_interpolated_transmissions)/airmass
    return extcoeff
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def interpolate_extcoeff(wave, wl_transm, extcoeff):
    interp_func = interp1d(wl_transm, extcoeff, kind='linear', fill_value='extrapolate')
    return interp_func(wave)
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def interpolate_eff(wave, wl, eff):
    interp_func = interp1d(wl, eff, kind='linear', fill_value='extrapolate')
    return interp_func(wave)
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def estimate_flux_correction_gaussian(seeing, fiber_diameter):
    # Full Width at Half Maximum (FWHM) of the Gaussian
    fwhm = seeing
    
    # Standard deviation of the Gaussian
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    
    # Radius of the fiber
    radius = fiber_diameter / 2
    
    # Flux within the fiber (analytical integration of the Gaussian from -radius to radius)
    flux_within_fiber = np.sqrt(np.pi) * sigma * (erf(radius / (np.sqrt(2) * sigma)) - erf(-radius / (np.sqrt(2) * sigma))) / 2

    # Total flux (scaling factor for a normalized Gaussian)
    total_flux = np.sqrt(np.pi) * sigma
    
    # Correction factor
    correction_factor =  flux_within_fiber / total_flux
    
    return correction_factor
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def interpolate_spectrum(wl_initial, spec_initial, wl_output):
    interpolator = interp1d(wl_initial, spec_initial, kind='linear', bounds_error=False, fill_value=0)
    spec_output = interpolator(wl_output)
    return spec_output
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def planck(wavelength, temperature):
    numerator = 2 * plank_constant * ligth_speed ** 2
    exponent = (plank_constant * ligth_speed) / (wavelength * boltzmann_constant * temperature)
    denominator = (wavelength ** 5) * (np.exp(exponent) - 1)
    black_body = numerator / denominator
    return black_body
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def safe_value(value, placeholder=-9999):
    if np.isnan(value) or np.isinf(value):
        return placeholder
    else:
        return value
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def mag_integral(y_values, x_values):
    if len(x_values) % 2 == 0:
        x_values = x_values[:-1]
        y_values = y_values[:-1]
    integral = simpson(y_values*x_values, x_values)
    #integral = simpson(y_values, x_values)
    #integral = np.sum(y_values)
    magnitude = -2.5*np.log10(integral)
    return magnitude
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def rolling_std(data, window_size):
    if window_size < 3:
        raise ValueError("Window size must be at least 3")
    pad_size = window_size // 2
    padded_data = np.pad(data, pad_size, mode='reflect')
    sq_data = padded_data**2
    mean = uniform_filter1d(padded_data, window_size)[pad_size:-pad_size]
    mean_sq = uniform_filter1d(sq_data, window_size)[pad_size:-pad_size]
    std = np.sqrt(np.clip(mean_sq - mean**2,0,None))
    return std
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def rebin_spectrum(wl, x, coarse_resolution):
    nnn = int(round((wl[-1] - wl[0]) / coarse_resolution))
    nwl = np.linspace(wl[0], wl[-1], nnn)
    spectrum_LR = np.zeros_like(nwl)

    for i in range(len(nwl)):
        start_idx = np.searchsorted(wl, nwl[i] - coarse_resolution / 2)
        end_idx = np.searchsorted(wl, nwl[i] + coarse_resolution / 2)
        spectrum_LR[i] = np.sum(x[start_idx:end_idx]) / (end_idx - start_idx)

    return nwl, spectrum_LR
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def flatten_spectrum(wl, x, coarse_resolution):
    nnn = int(round((wl[-1] - wl[0]) / coarse_resolution))
    nwl = np.linspace(wl[0], wl[-1], nnn)
    spectrum_LR = np.zeros_like(nwl)

    for i in range(len(nwl)):
        start_idx = np.searchsorted(wl, nwl[i] - coarse_resolution / 2)
        end_idx = np.searchsorted(wl, nwl[i] + coarse_resolution / 2)
        spectrum_LR[i] = np.sum(x[start_idx:end_idx]) / (end_idx - start_idx)

    f = interp1d(nwl, spectrum_LR, kind='linear')
    spectrum_envelope = f(wl)
    spectrum_flat_0 = x / spectrum_envelope
    spec_median = np.nanmedian(spectrum_flat_0)
    spectrum_flat = spectrum_flat_0 / spec_median

    rms = np.nanstd(spectrum_flat, ddof=1)
    skew = st.skew(spectrum_flat, nan_policy='omit')
    kurt = st.kurtosis(spectrum_flat, nan_policy='omit')

    return spectrum_flat, rms, skew, kurt
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def statistics(y):
    max_y = np.nanmax(y)
    min_y = np.nanmin(y)
    sum_y = np.nansum(y)
    mean_y = np.nanmean(y)
    median_y = np.nanmedian(y)
    rms_y = np.nanstd(y, ddof=1)
    skew_y = st.skew(y/max(y), nan_policy='omit')
    kurt_y = st.kurtosis(y/max(y), nan_policy='omit')
    return max_y, min_y, sum_y, mean_y, median_y, rms_y, skew_y, kurt_y
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
def bisector_analysis(wavelengths, flux, gk, sp, nps):
    smoothed_flux = gaussian_filter1d(flux, gk)

    npt=len(wavelengths)
    # If npt is not more than 1 or smoothed_flux is empty, return np.nan and 'U'
    if npt <= 1 or len(smoothed_flux) == 0:
        return [np.nan]*5 + ['U']

    max_flux = np.max(smoothed_flux)
    levels = np.linspace(max_flux, np.min(smoothed_flux), 100)
    bisector_wavelengths = []

    for level in levels:
        left_index = np.argmin(np.abs(smoothed_flux[:npt//2] - level))
        right_index = np.argmin(np.abs(smoothed_flux[npt//2:] - level)) + npt // 2
        bisector_wavelengths.append(0.5 * (wavelengths[left_index] + wavelengths[right_index]))
    
    bisector_wavelengths = np.array(bisector_wavelengths)
    
    # Compute statistics of bisector
    mean = np.mean(bisector_wavelengths)
    std = np.std(bisector_wavelengths)
    std = std if std > 0 else 1.0
    standardized_bisector = (bisector_wavelengths - mean) / std
    bisector_std = std
    
    bisector_skew = stats.skew(standardized_bisector)
    bisector_kurt = stats.kurtosis(standardized_bisector)
    center_wavelength = wavelengths[len(wavelengths) // 2]
    bisector_mean_position = np.mean(bisector_wavelengths)
    
    if bisector_mean_position < center_wavelength:
        bisector_asymmetry = "B"
    elif bisector_mean_position > center_wavelength:
        bisector_asymmetry = "R"
    else:
        bisector_asymmetry = "C"

    return bisector_wavelengths, bisector_mean_position, bisector_std, bisector_skew, bisector_kurt, bisector_asymmetry

# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def bisector_analysis1(wavelengths, flux, gk, sp, nps):
    smoothed_flux = gaussian_filter1d(flux, gk)

    # Performing spline interpolation
    #spline_flux = UnivariateSpline(wavelengths, smoothed_flux, k=sp, s=nps)
    #interpolation_function = interp1d(wavelengths, smoothed_flux, kind='cubic')
    #npt=3*len(wavelengths)
    #interpolated_wavelengths = np.linspace(wavelengths[0], wavelengths[-1], npt)
    #interpolated_flux = spline_flux(interpolated_wavelengths)
    #interpolated_flux = interpolation_function(interpolated_wavelengths)
    

    # Compute the bisector
    npt=len(wavelengths)
    max_flux = np.max(smoothed_flux)
    levels = np.linspace(max_flux, np.min(smoothed_flux), 100)
    bisector_wavelengths = []
    
    for level in levels:
        left_index = np.argmin(np.abs(smoothed_flux[:npt//2] - level))
        right_index = np.argmin(np.abs(smoothed_flux[npt//2:] - level)) + npt // 2
        bisector_wavelengths.append(0.5 * (wavelengths[left_index] + wavelengths[right_index]))
    
    bisector_wavelengths = np.array(bisector_wavelengths)
   
    
    # Compute statistics of bisector
    mean = np.mean(bisector_wavelengths)
    std = np.std(bisector_wavelengths)
    if std>0:
        True
    else:
        std=1.0
    standardized_bisector = (bisector_wavelengths - mean) / std
    bisector_std = std
    
    bisector_skew = stats.skew(standardized_bisector)
    bisector_kurt = stats.kurtosis(standardized_bisector)
    center_wavelength = wavelengths[len(wavelengths) // 2]
    bisector_mean_position = np.mean(bisector_wavelengths)
    
    if bisector_mean_position < center_wavelength:
        bisector_asymmetry="B"
        #print("The asymmetry is towards the blue (left) side.")
    elif bisector_mean_position > center_wavelength:
        bisector_asymmetry="R"
        #print("The asymmetry is towards the red (right) side.")
    else:
        bisector_asymmetry="C"
        #print("The line is symmetric with respect to its center.")
    #print(bisector_variance, bisector_skew, bisector_kurt, bisector_asymmetry)
    return bisector_wavelengths, bisector_mean_position, bisector_std, bisector_skew, bisector_kurt, bisector_asymmetry

# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def detect_peak_features(wavelengths, flux, height, gk, sp, nps):
    
    smoothed_flux = gaussian_filter1d(flux, gk)
    
    #window_length=5
    #polyorder=3
    #smoothed_flux = savgol_filter(flux, window_length, polyorder)

    
    # Create a spline interpolation of the flux
    #spline_flux = UnivariateSpline(wavelengths, smoothed_flux, k=sp, s=nps)
    
    # Compute the second derivative
    #second_derivative_flux = spline_flux.derivative(n=2)(wavelengths)
    #second_derivative_flux = smoothed_flux.derivative(n=2)(wavelengths)

    delta_wavelengths = wavelengths[1:] - wavelengths[:-1]
    delta_smoothed_flux = smoothed_flux[1:] - smoothed_flux[:-1]
    second_derivative_flux = (delta_smoothed_flux[1:] - delta_smoothed_flux[:-1]) / (delta_wavelengths[1:] + delta_wavelengths[:-1])
    
    # Make sure wavelengths array is consistent with second_derivative_flux array
    wavelengths_for_second_derivative = 0.5 * (wavelengths[2:] + wavelengths[:-2])



    # Count the number of positive and negative peaks in second derivative
    pos_peaks = 0
    neg_peaks = 0
    for i in range(1, len(second_derivative_flux) - 1):
        if second_derivative_flux[i-1] < second_derivative_flux[i] > second_derivative_flux[i+1]:
            pos_peaks += 1
        elif second_derivative_flux[i-1] > second_derivative_flux[i] < second_derivative_flux[i+1]:
            neg_peaks += 1

    # Classify based on the number of peaks
    peak_class = "UU"
    if height < 0:  # Emission
        if pos_peaks == 2:
            peak_class = "De" #"Double-peaked emission"
        elif pos_peaks == 1:
            peak_class = "Se" #"Single-peaked emission"
        elif pos_peaks == 0:
            peak_class = "Ue" #"Unclassified emission"
        elif pos_peaks >2:
            peak_class = "Me" #"Multi-peaked emission"
    else:  # Absorption
        if neg_peaks == 1:
            peak_class = "Ea" #"Absorption with emission core"
        elif neg_peaks == 2:
            peak_class = "Da" #"Absorption with emission core 2 peaks"
        elif neg_peaks > 2:
            peak_class = "Ma" #"Absorption with emission core >2 peaks"
        else:
            peak_class = "Ra" #"Regular absorption"
    return peak_class 
# ------------------------------  END   FUNCTION ------------------------------ #


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def fit_absorption_line(wavelengths, intensities, wc, ws, fit_lc):
    line_center0 = wc
    wl_span = ws
    mask = (wavelengths > (line_center0 - wl_span)) & (wavelengths < (line_center0 + wl_span))
    wl=wavelengths[mask]
    flux = intensities[mask]
    npt=len(wavelengths[mask])
    
    const_guess = np.mean(flux)
    amplitude_guess = np.min(flux) - const_guess
    mean_guess = line_center0
    stddev_guess = 0.1

    #fit Gaussian
    try:
        if fit_lc:
            popt_gaussian, pcov = curve_fit(negative_gaussian, wl, flux, p0=[const_guess, amplitude_guess, mean_guess, stddev_guess], maxfev=2000)
            y_fit_gaussian = negative_gaussian(wl, *popt_gaussian)
            errors = np.sqrt(np.diag(pcov))
            const_c = popt_gaussian[0]
            depth_c = popt_gaussian[1]
            lambda_c = popt_gaussian[2]
            sigma = np.abs(popt_gaussian[3])
            sigma_error = errors[3]
        else:
            popt_gaussian, pcov = curve_fit(lambda x, const, amplitude, stddev: negative_gaussian(x, const, amplitude, wc, stddev), wl, flux, p0=[const_guess, amplitude_guess, stddev_guess], maxfev=2000)
            y_fit_gaussian = negative_gaussian(wl, popt_gaussian[0], popt_gaussian[1], wc, popt_gaussian[2])
            errors = np.sqrt(np.diag(pcov))
            const_c = popt_gaussian[0]
            depth_c = popt_gaussian[1]
            lambda_c = wc
            sigma = np.abs(popt_gaussian[2])
            sigma_error = errors[2]

        ERROREG = 0
    except RuntimeError as e:
        ERROREG = 1
        sigma = np.nan
        sigma_error = np.nan

    # Fit Lorentzian
    ERRORE = 0
    try:
        # Get guess from gaussian fit
        if np.abs(const_c)>0 and np.abs(depth_c)>0 and np.abs(lambda_c)>0:
            const_guess = const_c
            amplitude_guess = depth_c
            mean_guess = lambda_c
            stddev_guess = sigma
        if fit_lc:
            popt_lorentzian, pcov = curve_fit(negative_lorentzian, wl, flux, p0=[const_guess, amplitude_guess, mean_guess, stddev_guess], maxfev=2000)
            y_fit_lorentzian = negative_lorentzian(wl, *popt_lorentzian)
            errors = np.sqrt(np.diag(pcov))  # Extract the errors from the covariance matrix
            continuum=popt_lorentzian[0]
            amplitude=popt_lorentzian[1]
            line_center=popt_lorentzian[2]
            gamma=np.abs(popt_lorentzian[3])
            continuum_error = errors[0]
            amplitude_error = errors[1]
            line_center_error = errors[2]
            gamma_error = errors[3]
        else:
            popt_lorentzian, pcov = curve_fit(lambda x, const, amplitude, stddev: negative_lorentzian(x, const, amplitude, wc, stddev), wl, flux, p0=[const_guess, amplitude_guess, stddev_guess], maxfev=2000)
            y_fit_lorentzian = negative_lorentzian(wl, popt_lorentzian[0], popt_lorentzian[1], wc, popt_lorentzian[2])
            errors = np.sqrt(np.diag(pcov))
            continuum = popt_lorentzian[0]
            amplitude = popt_lorentzian[1]
            line_center = wc
            gamma = np.abs(popt_lorentzian[2])
            continuum_error = errors[0]
            amplitude_error = errors[1]
            line_center_error = 0
            gamma_error = errors[2]
    except RuntimeError as e:
        ERRORE = 1
        continuum=np.nan
        amplitude=np.nan
        line_center=np.nan
        gamma=np.nan
        continuum_error=np.nan
        amplitude_error=np.nan
        line_center_error=np.nan
        gamma_error=np.nan
        reduced_chisq=np.nan
        r_squared=np.nan
        skewness=np.nan
        rms_rediduals=np.nan
        kurtosis=np.nan
        delta_wl=np.nan
        snr=np.nan
        rv=np.nan
        rv_error=np.nan
        area_line=np.nan
        area_lorentzian=np.nan
        blu_res=np.nan
        red_res=np.nan
        bisector=np.nan
        bisector_mean_position=np.nan
        bisector_variance=np.nan
        bisector_skew=np.nan
        bisector_kurt=np.nan
        bisector_std = np.nan
        bisector_asymmetry = "U"
        peak_features = "U"
        rms_residuals = np.nan


    if ERRORE == 0:
        delta_wl=line_center0-line_center
        rv = (line_center/line_center0 - 1)*const.c.to(u.km/u.s).value
        rv_error = const.c.to(u.km/u.s).value/line_center0 * line_center_error
        
        # Subset of flux and y_fit around ng*gamma
        ng = 3
        mask1 = (wl > (line_center - ng * gamma)) & (wl < (line_center + ng * gamma))
        wl_subset=wl[mask1]

        if wl_subset.size > 0:

            flux_subset = flux[mask1]/continuum
            y_fit_lorentzian_subset = y_fit_lorentzian[mask1]/continuum
            area_line       = np.trapz(flux_subset, x=wl_subset)
            area_lorentzian = np.trapz(y_fit_lorentzian_subset, x=wl_subset)
            # rms        
            residuals = flux_subset - y_fit_lorentzian_subset
            if residuals.size>0:
                rms_residuals = np.std(residuals)
            else:
                rms_residuals = np.nan
            if continuum >0:
                snr=abs(amplitude/continuum)/rms_residuals
            else:
                snr=np.nan
    
            # blue and red residuals
            # Find the index of the line center within the wl_subset
            center_index = np.argmin(np.abs(wl_subset - line_center))
            # Integrate the residuals for the left and right wings
            blu_res = np.trapz(residuals[:center_index + 1], x=wl_subset[:center_index + 1])/area_lorentzian
            red_res = np.trapz(residuals[center_index:], x=wl_subset[center_index:])/area_lorentzian
            #print(blu_res,red_res)
    
            # CHI2
            chisq = np.sum((residuals / np.std(flux_subset)) ** 2)
            dof = len(wl) - len(popt_lorentzian)
            reduced_chisq = chisq / dof
            # R2
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((flux_subset - np.mean(flux_subset)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)
                    
            flux1=1-flux_subset
            sum_flux = sum(flux1)
            npt_subset = len(flux_subset)
            skewness = sum(flux1[i] * (wl_subset[i] - line_center)**3 for i in range(npt_subset)) / (sum_flux * sigma**3)
            kurtosis = sum(flux1[i] * (wl_subset[i] - line_center)**4 for i in range(npt_subset)) / (sum_flux * sigma**4) - 3
    
            #new_wl = wl*(1-rv/const.c.to(u.km/u.s).value)
            #new_flux = flux
            #speck=new_flux/continuum
            #speck_wl=new_wl
            
            gk = 5 # number of points for gaussian kernel smoothing
            sp = 4 # polynomial degree for spline
            nps = 0 # The parameter nps helps to control the trade-off between closeness of fit and smoothness of the spline. 
            # If nps=0, the spline will interpolate through all data points (meaning it passes exactly through each point). 
            # Larger values of nps increase the amount of smoothing.
            
           
            bisector_results = bisector_analysis(wl_subset, flux_subset, gk, sp, nps)
            bisector_mean_position=bisector_results[1]
            bisector_std=bisector_results[2]
            bisector_skew=bisector_results[3]
            bisector_kurt=bisector_results[4]
            bisector_asymmetry=bisector_results[5]
    
            peak_features = detect_peak_features(wl_subset, flux_subset, amplitude, gk, sp, nps)

        else:
            area_line = np.nan
            area_lorentzian = np.nan
            skewness = np.nan
            kurtosis = np.nan
            bisector_mean_position = np.nan
            bisector_std = np.nan
            bisector_skew = np.nan
            bisector_kurt = np.nan
            bisector_asymmetry = "U"
            peak_features = "U"
            rms_residuals = np.nan
            snr = np.nan
            reduced_chisq = np.nan
            r_squared = np.nan
            blu_res=np.nan
            red_res=np.nan
        
    else:
        nothing = True

    results=[continuum,continuum_error,amplitude,amplitude_error,line_center,line_center_error,gamma,gamma_error,sigma,sigma_error,
    skewness,kurtosis,rms_residuals,snr,area_line,area_lorentzian,blu_res,red_res,
    reduced_chisq,r_squared,delta_wl,rv,rv_error,bisector_mean_position,bisector_std,bisector_skew,bisector_kurt,bisector_asymmetry,peak_features]
    
    
    return results
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

    # no atmospheric refraction
    star_altaz = star_coords.transform_to(altaz_frame)
    star_azimuth0 = star_altaz.az
    star_altitude0 = star_altaz.alt

    # with atmospheric refraction
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

    sun_ra=sun_ra_dec.ra.deg
    sun_dec=sun_ra_dec.dec.deg
    sun_altitude = sun_altaz.alt.deg
    sun_azimuth = sun_altaz.az.deg

    return sun_ra, sun_dec, sun_altitude, sun_azimuth
# ------------------------------  END   FUNCTION ------------------------------ #

# ------------------------------  BEGIN FUNCTION ------------------------------ #
def calculate_moon_position(mjd, lon, lat, elevation):
    location = EarthLocation(lon=lon*u.deg, lat=lat*u.deg, height=elevation*u.m)
    time = Time(mjd, format='mjd')

    moon_altaz = get_moon(time).transform_to(AltAz(obstime=time, location=location))
    moon_ra_dec = get_moon(time).transform_to('icrs')

    moon_ra=moon_ra_dec.ra.deg
    moon_dec=moon_ra_dec.dec.deg
    moon_altitude = moon_altaz.alt.value
    moon_azimuth = moon_altaz.az.value

    return moon_ra, moon_dec, moon_altitude, moon_azimuth
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


# ------------------------------  BEGIN FUNCTION ------------------------------ #
def convert_arabic_to_uppercase_roman_numbers(arabic):
    roman_numerals = {1: 'I', 4: 'IV', 5: 'V', 9: 'IX', 10: 'X', 40: 'XL', 50: 'L', 90: 'XC', 100: 'C', 400: 'CD', 500: 'D', 900: 'CM', 1000: 'M'}
    roman = ''
    for value, numeral in sorted(roman_numerals.items(), reverse=True):
        while arabic >= value:
            arabic -= value
            roman += numeral
    return roman.upper()
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


variables_fit_line = ['continuum','continuum_error','amplitude','amplitude_error','line_center','line_center_error','gamma','gamma_error','sigma','sigma_error','skewness','kurtosis','rms_residuals','snr','area_line','area_lorentzian','blu_res','red_res','reduced_chi2','r2','delta_wl','rv','rv_error','bisector_mean_position','bisector_variance','bisector_skew','bisector_kurt','bisector_asymmetry','peak_features']

# CCD wavelengths
wl_min_b = 3783.0
wl_max_b = 5304.0
wl_min_r = 5338.0
wl_max_r = 6912.0


# ==============================  BEGIN BLOCK  ================================ #
hdu = fits.open('phoenix_halpha_spectrum_slices.fits')
pha = hdu[1].data
hdu.close()
wl_pha = np.array(pha['spectrum_slice'][0][:])
pha_nspec = len(pha)

pha_meh_criteria            = np.isin(pha['meh'], [-2, -1, 0])
pha_logg_criteria           = np.isin(pha['logg'], [0, 1, 2, 3, 3.5, 4, 4.5, 5, 6])
pha_logte_abs_criteria1     = np.any([np.abs(10 ** pha['logte'] - i) < 50 for i in range(3000, 10001, 1000)], axis=0)
pha_logte_abs_criteria2     = np.any([np.abs(10 ** pha['logte'] - i) < 50 for i in range(3500, 7501, 1000)], axis=0)
pha_logte_gte_3000_criteria = 10 ** pha['logte'] >= 3000

# Combine the criteria
pha_combined_criteria = pha_meh_criteria & pha_logg_criteria & (pha_logte_abs_criteria1 | pha_logte_abs_criteria2) & pha_logte_gte_3000_criteria
# Apply the selection criteria
pha_selected = pha[pha_combined_criteria]

#pha_selected_indices = np.arange(pha_nspec)[pha_combined_criteria]
# Collect the spectra for the selected indices
#pha_selected_spectra = np.array([np.array(pha['spectrum_slice'][i][:]) for i in pha_selected_indices])
# ==============================  END   BLOCK  ================================ #


# ==============================  BEGIN BLOCK  ================================ #
# Lines
with fits.open('spectral_line_list_coluzzi_harps.fits') as absorption_lines:
    absorption_lines_data_coluzzi = absorption_lines[1].data

mask_common_lines = np.array([len(comment) > 0 for comment in absorption_lines_data_coluzzi['COMMENT']])
#mask_common_lines = np.array([(len(comment) > 0 or len(sptypes) > 0) for comment, sptypes in zip(absorption_lines_data_coluzzi['COMMENT'], absorption_lines_data_coluzzi['SpTypes'])])
common_lines = absorption_lines_data_coluzzi[mask_common_lines]
# ==============================  END   BLOCK  ================================ #



# ==============================  BEGIN BLOCK  ================================ #
with fits.open('nist_lines_output.fits') as absorption_lines:
    absorption_lines_data = absorption_lines[1].data
# ==============================  END   BLOCK  ================================ #


# ==============================  BEGIN BLOCK  ================================ #
nist_selected_lines = pd.read_csv('nist_selected_lines.csv')
# ==============================  END   BLOCK  ================================ #


# ==============================  BEGIN BLOCK  ================================ #
wl_main_lines = [
    (3933.664, 10, 'CaH'),
    (3970.074, 10, 'Hepsilon'),
    (4101.737, 10, 'Hdelta'),
    (4340.468, 10, 'Hgamma'),
    (4471.477, 10, 'HeI'),
    (4481.129, 10, 'MgII'),
    (4541.590, 10, 'HeII'),
    (4685.682, 10, 'HeII'),
    (4861.332, 10, 'Hbeta'),
    (5270.39,  10, 'Fe'),
    (5895.92,  10, 'NaI_D1'),
    (5889.95,  10, 'NaI_D2'),
    (6707.815, 10, 'LiI')
]
# ==============================  END   BLOCK  ================================ #

# ==============================  BEGIN BLOCK  ================================ #
wl_important_lines = [
    (4861.332, 10, 'Hbeta'),
    (4340.468, 10, 'Hgamma'),
    (4101.737, 10, 'Hdelta'),
    (3970.074, 10, 'Hepsilon'),
    (4471.477, 10, 'HeI'),
    (4541.590, 10, 'HeII'),
    (3933.664, 10, 'CaH'),
    (4481.129, 10, 'MgII'),
    (5269.541, 10, 'FeI'),
    (5895.923, 10, 'NaI_D1'),
    (5889.953, 10, 'NaI_D2'),
    (6707.815, 10, 'LiI'),
    (4554.033, 10, 'BaII'),
    (4262.270, 10, 'TcI')
]
# ==============================  END   BLOCK  ================================ #

# ==============================  BEGIN BLOCK  ================================ #
print_filters = False
if print_filters == True:
    #    (3968.470, 'CaK'),
    wavelengths = [
        (3797.900, 'H10'),
        (3835.386, 'H9'),
        (3889.051, 'H8'),
        (3933.664, 'CaH'),
        (3970.074, 'Heps'),
        (4101.737, 'Hdel'),
        (4340.468, 'Hgam'),
        (4471.477, 'HeI'),
        (4481.129, 'MgII'),
        (4541.590, 'HeII'),
        (4685.682, 'HeII'),
        (4861.332, 'Hbet'),
        (5892.938, 'NaI'),
        (6562.817, 'Half'),
        (6707.815,' LiI')
    ]
    
    for i, (wavelength, label) in enumerate(wavelengths):
        lower_bound = max(wl_min_b, wavelength - 10)
        upper_bound = min(wl_max_r, wavelength + 10)
        print(f"{lower_bound:8.3f} {upper_bound:8.3f} | {wavelength:8.3f} {label}")
    
        if i < len(wavelengths) - 1:
            next_wavelength, next_label = wavelengths[i + 1]
            
            interval_start = min(wl_max_r, wavelength + 10)
            interval_end = max(wl_min_b, next_wavelength - 10)
            
            # Check if the interval contains wl_max_b and wl_min_r
            if interval_start < wl_max_b and interval_end > wl_max_b:
                print(f"{interval_start:8.3f} {wl_max_b:8.3f} | {wl_max_b - wavelength:8.3f}")
                print(f"{wl_min_r:8.3f} {interval_end:8.3f} | {next_wavelength - wl_min_r:8.3f}")
            else:
                print(f"{interval_start:8.3f} {interval_end:8.3f} | {interval_end - interval_start:8.3f}")
    
    # Extra line after Half
    if wavelengths[-1][1] == 'Half':
        print(f"{wavelengths[-1][0] + 10:8.3f} {wl_max_r:8.3f} | {wl_max_r - wavelengths[-1][0]:8.3f}")
    
    print()
    print()

    w=5
    for i, (wavelength, label) in enumerate(wavelengths):
        print(f"{wavelength - w:8.3f} {wavelength + w:8.3f} | {wavelength:8.3f} {label}")
        if i < len(wavelengths) - 1:
            next_wavelength, next_label = wavelengths[i + 1]
            print(f"{wavelength + w:8.3f} {next_wavelength - w:8.3f} | {next_wavelength - wavelength:8.3f}")
    
    print()
    print()
    
    # Minimum and maximum usable wavelengths in the blue and red CCDs in Angstrom
    wl_min_b = 3783.0
    wl_max_b = 5304.0
    bw_b = wl_max_b - wl_min_b  # 1521.0
    wl_min_r = 5338.0
    wl_max_r = 6912.0
    bw_r = wl_max_r - wl_min_r  # 1574.0
    
    
    num_filters = 2
    filter_width_b = (wl_max_b - wl_min_b) / num_filters
    filter_width_r = (wl_max_r - wl_min_r) / num_filters
    
    # Iterate over filters in the blue range
    for i in range(num_filters):
        filter_min = wl_min_b + i * filter_width_b
        filter_max = filter_min + filter_width_b
        print("BLU", i+1, filter_min, filter_max)
        
        # Iterate through rows and print desired columns if lambda is in range
        for line in common_lines:
            if filter_min <= line['lambda'] <= filter_max:
                print(line['lambda'], line['element'], line['COMMENT'], line['SpTypes'])
        print("\n")
    
    # Iterate over filters in the red range
    for i in range(num_filters):
        filter_min = wl_min_r + i * filter_width_r
        filter_max = filter_min + filter_width_r
        print("RED",i+1,filter_min, filter_max)
        
        # Iterate through rows and print desired columns if lambda is in range
        for line in common_lines:
            if filter_min <= line['lambda'] <= filter_max:
                print(line['lambda'], line['element'], line['COMMENT'], line['SpTypes'])
        print("\n")
    sys.exit()
# ==============================  END   BLOCK  ================================ #



# ==============================  BEGIN BLOCK  ================================ #
filename = 'filters.fits'
hdul = fits.open(filename)
filter_names = [
    "TYCHO.B_bes",
    "TYCHO.V_bes",
    "Hipparcos.Hp_bes",
    "SDSS.g",
    "SDSS.r",
    "GAIA3.Gbp"
]
filters_data = {}
for filter_name in filter_names:
    ext_index = None
    for i, hdu in enumerate(hdul):
        if hdu.name == filter_name:
            ext_index = i
            break
    if ext_index is not None:
        data = hdul[ext_index].data
        wavelength = data.field(0)
        transmission = data.field(1)
        max_transmission = np.max(transmission)
        normalized_transmission = transmission / max_transmission
        filters_data[filter_name] = {'wavelength': wavelength, 'transmission': normalized_transmission}
hdul.close()
# ==============================  END   BLOCK  ================================ #


# ==============================  BEGIN BLOCK  ================================ #
trasm_df = pd.read_csv('smoothed_atmospheric_transmissions.csv')
wl_trasm = trasm_df['wavelength']
# ==============================  END   BLOCK  ================================ #


# ==============================  BEGIN BLOCK  ================================ #
eff_df = pd.read_csv('smoothed_instrument_eff.csv')
wl_eff = eff_df['wavelength']
total_eff=eff_df['telescope']*eff_df['fiber']*eff_df['instrument']*eff_df['ccd']
# ==============================  END   BLOCK  ================================ #



# ==============================  BEGIN BLOCK  ================================ #
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
# ==============================  END   BLOCK  ================================ #

    


# ==============================  BEGIN BLOCK  ================================ #
file_names = ['adp', 'ccf', 'bis', 'gui']
base_filename = '{}_small.hdr'
spectral_masks = ["G2", "K0", "K5", "M2", "M4"]
# Dictionary to store keywords and comments for different file types
datakw = {}
# Loop through each file type containing the keywords
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
# ==============================  END   BLOCK  ================================ #


# ==============================  BEGIN BLOCK  ================================ #
# READ HARPS CHARACTERISTICS
hdu_harps = fits.open('harps_output.fits', ignore_missing_end=True)
harps_data = hdu_harps[1].data
hdu_harps.close()
harps_wavelength0=harps_data['wavelength']
harps_trasm_atm0 =harps_data['throughput_atmosphere']
harps_trasm_det0 =harps_data['throughput_detector']
harps_trasm_fib0 =harps_data['throughput_fibinj']
harps_trasm_inst0=harps_data['throughput_instrument']
harps_trasm_tel0 =harps_data['throughput_telescope']
harps_trasm_bla0 =harps_data['throughput_blaze']
harps_dispersion0=harps_data['dispersion_dispersion']
harps_sed_sky0   =harps_data['sed_sky']
harps_trasm_all  =harps_trasm_det0*harps_trasm_fib0*harps_trasm_inst0*harps_trasm_tel0
#
# ==============================  END   BLOCK  ================================ #


# ==============================  BEGIN BLOCK  ================================ #
# READ TABLE WITH SEEING VALUES AS FUNCTION OF TIME
filename = 'table_seeing.fits'
hdul = fits.open(filename)
tabseeing = hdul[1].data
mjd_seeing=tabseeing['mjd']
tab_seeing=tabseeing['fwhm']
hdul.close()
# ==============================  END   BLOCK  ================================ #



#********************************************************************************
#
#  M A I N   C O D E
#
#********************************************************************************




# Open the catalog FITS file
with fits.open('harps_simbad_crossmatch_nearest.fits') as harps_list_obs:
    harps_list_obs_data = harps_list_obs[1].data
    indices = np.arange(len(harps_list_obs_data))
    # Shuffle the indices
    np.random.shuffle(indices)
    # Use the shuffled indices to access the rows in random order
    shuffled_harps_list_obs_data = harps_list_obs_data[indices]

masking = False
if masking == True:
    print('begin masking')
    
    mask_spt = [not row['spt_obj'].startswith(('O','B','A','F', 'G', 'K', 'M', 'N', '-')) for row in shuffled_harps_list_obs_data]
    #mask_spt = [row['spt_obj'].startswith(('C','W','Q','R', 'S', '?')) for row in shuffled_harps_list_obs_data]
    #mask_acen = [(row['ra']>219.84) & (row['ra']<219.9) & (row['dec']<-60.826) & (row['dec']>-60.844)  for row in harps_list_obs_data]
    #mask_bethyi= [(row['ra']>6.34) & (row['ra']<6.48) & (row['dec']>-77.265) & (row['dec']<-77.245)  for row in harps_list_obs_data]
    #mask_pole= [(row['dec']>-70) & (row['dec']<-90)  for row in harps_list_obs_data]
    
    #mask=mask_acen
    #mask=mask_bethyi
    #mask=mask_pole
    mask=mask_spt
    
    ## Apply the mask to select only the rows that meet the conditions
    #selected_rows = [row for row, m in zip(harps_list_obs_data, mask) if m]
    selected_rows = [row for row, m in zip(shuffled_harps_list_obs_data, mask) if m]
    
    ## Sort the selected_rows by 'mjd_obs' in descending order
    selected_rows = sorted(selected_rows, key=lambda x: x['mjd_obs'], reverse=False)
    
    #print(len(selected_rows))
    #for row in selected_rows:
    #    print(row['spt_obj'])
        
    print('end masking')

#sys.exit()

k = 0

#speck = None
#speck_wl = None



#sys.exit()


# prepare the csv output
# Read column names from the ASCII file
with open('harps_extract_all_data.output.hdr.csv', 'r') as out_hdr_file:
    out_header_line = out_hdr_file.readline().strip()
    out_column_names = out_header_line.split(',')

outfile = 'harps_extract_all_data.output.csv'

with open(outfile, 'w', newline='\n') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(out_column_names)
#

start_time = datetime.datetime.now()
previous_time = start_time

# all
numtest=323653

#numtest=2
# Loop over all rows


# here uncomment the correct for loop according to the test to be performed


#for row in harps_list_obs_data:
### loop on all the data as stored in the table

#for row in selected_rows:
### loop on a sample of selected data

for row in shuffled_harps_list_obs_data:
## loop on random sorted data

    k += 1
    #print(k)
    if(k>numtest+1):
        sys.exit()
    dp_id_raw = row['dp_id_raw']
    archive_id_spectra0 = row['archive_id_spectra'].strip()
    archive_id_ancillary0 = row['archive_id_ancillary'].strip()

    if len(archive_id_spectra0) > 0:
        archive_id_spectra = archive_id_spectra0 + '.fits'
        dateobs = dp_id_raw[6:16]
        spectra_file = os.path.join(storage_path, dateobs, archive_id_spectra)
        all_results = {}

# ++++++++++++++++++++++++  BEGIN ADP SPECTRA  +++++++++++++++++++++++ #
        adpisfits = testfilefits(spectra_file)
        filenames_dict = defaultdict(dict)
        filenames_dict['RAW']['value'] = dp_id_raw
        filenames_dict['ADP']['value'] = archive_id_spectra0
        filenames_dict['ADP']['status'] = adpisfits

        if adpisfits == True:
            harps_obs = fits.open(spectra_file)
            current_time = datetime.datetime.now()
            elapsed_total = current_time - start_time
            elapsed_previous = current_time - previous_time
            #print(f'adp {spectra_file.replace("/scratch/mbarbier/mylocal/machines/Desktop_media/SeagateHub/DATA/HARPS/harpstar/data/","")}')
            print(f'{k}, {spectra_file.replace("/scratch/mbarbier/mylocal/machines/Desktop_media/SeagateHub/DATA/HARPS/harpstar/data/","")}, {current_time}, {elapsed_total}, {elapsed_previous}')
            previous_time = current_time
            hdr0 = harps_obs[0].header
            hdr1 = harps_obs[1].header
            spec = harps_obs[1].data
          

            adp_data = {}
            keywords = datakw['adp']['keywords']
            for keyword in keywords:
                keyword_value = hdr0.get(keyword, None)
                #keyword_comment = hdr0.comments[keyword]
                try:
                    keyword_comment = hdr0.comments[keyword]
                except KeyError:
                    keyword_comment = ''
                column_name = keyword.replace('HIERARCH ESO ', '').replace(' ', '_').lower()
                adp_data[column_name] = {'value': keyword_value,'comment': keyword_comment}
                #print(adp_data[column_name])
            
            
            adp_dati = defaultdict(dict)
            adp_dati['idx']['value'] = k
            adp_dati['p3file']['value'] = spectra_file

            # read spectrum data
            wl = spec.field('WAVE')
            flux = spec.field('FLUX')
            npoints = flux.size
            # convert to 1D vectors
            wl_ = spec.field('WAVE').flatten()
            flux_ = spec.field('FLUX').flatten()

            #adp_data['tel_targ_delta_epoch'] = {'value': None,'comment': None}
            #adp_data['tel_targ_alpha_epoch'] = {'value': None,'comment': None}
            #adp_data['mjd-obs_center'] = {'value': None,'comment': None}
            #adp_data['diff_ra_ep_pnt'] = {'value': None,'comment': None}
            #adp_data['diff_de_ep_pnt'] = {'value': None,'comment': None}
            #adp_data['diff_ra_ob_pnt'] = {'value': None,'comment': None}
            #adp_data['diff_de_ob_pnt'] = {'value': None,'comment': None}
            #adp_data['num_pmra'] = {'value': None,'comment': None}
            #adp_data['num_pmde'] = {'value': None,'comment': None}
            #adp_data[''] = {'value': None,'comment': None}

            #j2000=51544.5
            coord_epoch=adp_data['tel_targ_epoch']['value']
            coord_epoch_mjd=Time(coord_epoch, format='decimalyear').mjd
#            if coord_epoch<1991:
#                print("*****  STRANGE EPOCH !!!!!    ",coord_epoch)
            mjd=adp_data['mjd-obs']['value']
            tmmean = adp_data['ins_det1_tmmean']['value']
            adp_dati['mjd_obs_center']['value']=mjd+adp_data['exptime']['value']*tmmean
            year_obs=Time(mjd, format='mjd').decimalyear
            adp_dati['year_obs']['value'] = year_obs

            tel_targ_alpha_degrees, tel_targ_delta_degrees, tel_targ_alpha_sexag, tel_targ_delta_sexag = convert_from_pseudosexagesimal_to_degrees(adp_data['tel_targ_alpha']['value'],adp_data['tel_targ_delta']['value'])
            adp_data['tel_targ_alpha']['value'] = tel_targ_alpha_degrees
            adp_data['tel_targ_delta']['value'] = tel_targ_delta_degrees

            ra1 , de1, _ , _ = convert_from_pseudosexagesimal_to_degrees(adp_data['ins_adc1_ra']['value'],adp_data['ins_adc1_dec']['value'])
            adp_data['ins_adc1_ra']['value']  = ra1
            adp_data['ins_adc1_dec']['value'] = de1

            ra2 , de2, _ , _ = convert_from_pseudosexagesimal_to_degrees(adp_data['ins_adc2_ra']['value'],adp_data['ins_adc2_dec']['value'])
            adp_data['ins_adc2_ra']['value']  = ra2
            adp_data['ins_adc2_dec']['value'] = de2
            
            pmra=adp_data['tel_targ_pma']['value']
            pmde=adp_data['tel_targ_pmd']['value']
            tel_targ_alpha_epoch = tel_targ_alpha_degrees + pmra/3600*(mjd-coord_epoch_mjd)/365.25
            tel_targ_delta_epoch = tel_targ_delta_degrees + pmde/3600*(mjd-coord_epoch_mjd)/365.25
            adp_dati['tel_targ_alpha_epoch']['value'] = tel_targ_alpha_epoch
            adp_dati['tel_targ_delta_epoch']['value'] = tel_targ_delta_epoch
            pm_tot_ob = np.sqrt(pmra**2+pmde**2)
            adp_dati['pm_tot_ob']['value']  = pm_tot_ob

            diff_ra_ep_pnt = (adp_dati['tel_targ_alpha_epoch']['value']-adp_data['ra']['value'])*3600*np.cos(adp_data['tel_targ_delta']['value']*np.pi/180)
            diff_de_ep_pnt = (adp_dati['tel_targ_delta_epoch']['value']-adp_data['dec']['value'])*3600
            diff_ra_ob_pnt = (adp_data['tel_targ_alpha']['value']-adp_data['ra']['value'])*3600*np.cos(adp_data['tel_targ_delta']['value']*np.pi/180)
            diff_de_ob_pnt = (adp_data['tel_targ_delta']['value']-adp_data['dec']['value'])*3600
            adp_dati['diff_ra_ep_pnt']['value']  = diff_ra_ep_pnt
            adp_dati['diff_de_ep_pnt']['value']  = diff_de_ep_pnt
            adp_dati['diff_ra_ob_pnt']['value']  = diff_ra_ob_pnt
            adp_dati['diff_ra_ob_pnt']['value']  = diff_de_ob_pnt

            #altazimutal coordinates
            azimuth_angle = true_azimuth(adp_data['tel_az']['value'])
            adp_data['tel_az']['value'] = azimuth_angle
            altitude_angle=adp_data['tel_alt']['value']
            
            # telescope coordinates
            #https://www.eso.org/sci/facilities/lasilla/telescopes/3p6/overview.html
            #The telescope is at a geographical location of:
            #70° 43' 54.1" W  -29° 15' 39.5" S (WGS84)
            #2400 metres above sea level
            lon=-70.731694
            lat=-29.260972
            ele=2400.0
            adp_data['tel_geolon']['value']=lon
            adp_data['tel_geolat']['value']=lat
            adp_data['tel_geoelev']['value']=ele
            
            # atmospheric parameters for refraction calculation
            humidity=adp_data['tel_ambi_rhum']['value']
            temperature=adp_data['tel_ambi_temp']['value']
            pressure=adp_data['tel_ambi_pres_start']['value']

            # sun/moon coordinates and distances from target
            sun_ra , sun_dec , sun_alt , sun_az  = calculate_sun_position(mjd, lon, lat, ele)
            moon_ra, moon_dec, moon_alt, moon_az = calculate_moon_position(mjd, lon, lat, ele)
            sun_coord  = SkyCoord(ra=sun_ra * u.deg, dec=sun_dec * u.deg, frame='icrs')
            moon_coord = SkyCoord(ra=moon_ra * u.deg, dec=moon_dec * u.deg, frame='icrs')
            object_coord = SkyCoord(ra=tel_targ_alpha_epoch*u.deg, dec=tel_targ_delta_epoch*u.deg, frame='icrs')
            angular_distance_obj_sun = object_coord.separation(sun_coord)
            angular_distance_obj_moon = object_coord.separation(moon_coord)
            #print(angular_distance_obj_sun.deg,angular_distance_obj_moon.deg,sun_alt.deg)
            adp_dati['sun_ra']['value']  =     sun_ra
            adp_dati['sun_dec']['value']  =    sun_dec
            adp_dati['sun_alt']['value']  =    sun_alt
            adp_dati['sun_az']['value']  =     sun_az
            adp_dati['moon_ra']['value']  =    moon_ra
            adp_dati['moon_dec']['value']  =   moon_dec
            adp_dati['moon_alt']['value']  =   moon_alt
            adp_dati['moon_az']['value']  =    moon_az
            adp_dati['angular_distance_obj_sun']['value']  =  angular_distance_obj_sun.deg
            adp_dati['angular_distance_obj_moon']['value']  = angular_distance_obj_moon.deg

            # calculate the azimuth and altitude using a full refraction model
            azimuth, altitude, azimuth_no_atm, altitude_no_atm = calculate_star_altaz(tel_targ_alpha_epoch, tel_targ_delta_epoch, lon, lat, ele, mjd, humidity, temperature, pressure)
            diff_az  = (azimuth.deg-adp_data['tel_az']['value'])*3600*np.cos(altitude*np.pi/180)
            diff_alt = (altitude.deg-adp_data['tel_alt']['value'])*3600
            adp_dati['azimuth']['value']  = azimuth.deg
            adp_dati['altitude']['value']  = altitude.deg
            adp_dati['azimuth_no_atm']['value']  = azimuth_no_atm.deg
            adp_dati['altitude_no_atm ']['value']  = altitude_no_atm.deg 
            adp_dati['diff_az']['value']  = diff_alt
            adp_dati['diff_alt']['value']  = diff_az

            #misc parameters
            #snr = adp_data['snr']['value']
            #cttot = adp_data['ins_det1_cttot']['value']
            exptime = adp_data['exptime']['value']
            fluxtot=np.sum(flux_)
            adp_dati['flux_total']['value']  = fluxtot

            # Spectral Type
            dpr_value = adp_data['dpr_type']['value']
            split_values = dpr_value.split(',')
            if len(split_values) >= 3:
                spt = split_values[2]
            else:
                spt = "NONE"
            adp_dati['spt_dpr']['value']  = spt

            # airmass
            airmass = (adp_data['tel_airm_start']['value']+adp_data['tel_airm_end']['value'])/2
            #airmass = adp_data['tel_airm_start']['value']
            if airmass>2.9:
                airmass_=2.9
            if airmass<1:
                airmass_=1
            adp_dati['airmass_start']['value']  = airmass

            #extinction coefficients 
            extcoeff = calc_extcoeff(airmass)
            extcoeff_interp = interpolate_extcoeff(wl_, wl_trasm, extcoeff)
            eff_interp = interpolate_eff(wl_, wl_eff, total_eff)

            # Seeing from header
            seeing_end=adp_data['tel_ambi_fwhm_end']['value']
            seeing_start=adp_data['tel_ambi_fwhm_start']['value']
            # Create a mask for the MJD values within the specified range
            mask_within_range = (mjd_seeing >= mjd) & (mjd_seeing <= mjd + exptime / 86400) & (tab_seeing>0)
            # Extract the seeing values within the specified MJD range
            seeing_within_range = tab_seeing[mask_within_range]
            mjd_seeing_within_range = mjd_seeing[mask_within_range]
            # Check if there is at least two elements within range
            ssize = seeing_within_range.size
            seeing_npt=ssize
            seeing_trend = 0
            if ssize > 2 :
                #print('seeing:', seeing_within_range)
                seeing_mean = np.mean(seeing_within_range)
                seeing_std = np.std(seeing_within_range)
                try:
                    coeff_fit_seeing, cov_fit_seeing = np.polyfit(mjd_seeing_within_range, seeing_within_range, 1, cov=True)
                    slope_seeing, intercept_seeing = coeff_fit_seeing
                    slope_error_seeing, intercept_error_seeing = np.sqrt(np.diag(cov_fit_seeing))
                    if np.abs(slope_seeing) > slope_error_seeing:
                        if slope_seeing > 0:
                            trend_seeing = -1 #-1  # degrading
                        elif slope_seeing < 0:
                            trend_seeing = 1  # improving
                except np.linalg.LinAlgError:
                    seeing_trend = 0
            elif ssize == 2 :
                seeing_mean=np.mean(seeing_within_range)
                seeing_std=0.2
                if seeing_within_range[1]>seeing_within_range[0]:
                    seeing_trend = -1
                elif seeing_within_range[1]<seeing_within_range[0]:
                    seeing_trend = 1
                else:
                    seeing_trend = 0
            elif ssize == 1 :
                seeing_mean=seeing_within_range[0]
                seeing_std=0.2
                seeing_trend = 0
            else:                    
                if(seeing_end>0 and seeing_start>0):
                    seeing_npt=2
                    seeing_mean=(seeing_start+seeing_end)/2
                    seeing_std=np.abs(seeing_start-seeing_end)/2
                    if seeing_end>seeing_start:
                        seeing_trend = -1
                    elif seeing_end<seeing_start:
                        seeing_trend = 1
                    else:
                        seeing_trend = 0
                elif(seeing_end<0 and seeing_start>0):
                    seeing_npt=1
                    seeing_mean=seeing_start
                    seeing_std=0.2
                    seeing_trend = 0
                elif(seeing_end>0 and seeing_start<0):
                    seeing_npt=1
                    seeing_mean=seeing_end
                    seeing_std=0.2
                    seeing_trend = 0
                else:
                    seeing_npt=0
                    seeing_mean=-1
                    seeing_std=-1
                    seeing_trend = 0
                    
            seeing_flag=1
            if seeing_mean<0:
                seeing_flag = -1
            adp_dati['seeing_mean']['value']  = seeing_mean
            adp_dati['seeing_std']['value']   = seeing_std
            adp_dati['seeing_npt']['value']   = seeing_npt
            adp_dati['seeing_trend']['value'] = seeing_trend
            adp_dati['seeing_flag']['value']  = seeing_flag
            

            # size of the atmospheric dispersion
            atmdisp = calc_atm_dispersion(altitude,temperature,pressure,humidity,year_obs)
            
            if seeing_flag == True:
                composite_fwhm = np.sqrt(seeing_mean**2 + atmdisp**2)
            else:
                composite_fwhm = atmdisp
            adp_dati['composite_fwhm']['value']  = composite_fwhm
            adp_dati['atmospheric_dispersion']['value']  = atmdisp
                        
            fiber_diameter = 1.0 # diameter of the fiber in arcseconds
            # calculate how many flux is loss due to seeing and atmospheric dispersion
            flux_seeing_correction=estimate_flux_correction_gaussian(composite_fwhm, fiber_diameter)
            adp_dati['flux_seeing_correction']['value']  = flux_seeing_correction
            # correct the flux           
            flux_corr=flux_/(exptime*flux_seeing_correction)*(extcoeff_interp/eff_interp)
            #calculate the magnitudes
            magnitudes,fluxes = clean_flux_and_calculate_magnitudes(flux_corr, wl_, filters_data)
            hpmag = magnitudes.get("Hipparcos.Hp_bes")
            btmag = magnitudes.get("TYCHO.B_bes")
            vtmag = magnitudes.get("TYCHO.V_bes")
            gmag  = magnitudes.get("SDSS.g")
            rmag  = magnitudes.get("SDSS.r")
            bpmag = magnitudes.get("GAIA3.Gbp")
            
            hpflux = fluxes.get("Hipparcos.Hp_bes")
            btflux = fluxes.get("TYCHO.B_bes")
            vtflux = fluxes.get("TYCHO.V_bes")
            gflux  = fluxes.get("SDSS.g")
            rflux  = fluxes.get("SDSS.r")
            bpflux = fluxes.get("GAIA3.Gbp")

            adp_dati['hpmag']['value']   = hpmag
            adp_dati['btmag']['value']   = btmag
            adp_dati['vtmag']['value']   = vtmag
            adp_dati['gmag']['value']    = gmag 
            adp_dati['rmag']['value']    = rmag 
            adp_dati['bpmag']['value']   = bpmag
            adp_dati['hpflux']['value']  = hpflux
            adp_dati['btflux']['value']  = btflux
            adp_dati['vtflux']['value']  = vtflux
            adp_dati['gflux']['value']   = gflux 
            adp_dati['rflux']['value']   = rflux 
            adp_dati['bpflux']['value']  = bpflux
            
            # don't apply the flux correction and calculates the magnitudes
            flux_corr=flux_/(exptime)*(extcoeff_interp/eff_interp)
            magnitudes,fluxes = clean_flux_and_calculate_magnitudes(flux_corr, wl_, filters_data)

            hpmag_nc = magnitudes.get("Hipparcos.Hp_bes")
            btmag_nc = magnitudes.get("TYCHO.B_bes")
            vtmag_nc = magnitudes.get("TYCHO.V_bes")
            gmag_nc  = magnitudes.get("SDSS.g")
            rmag_nc  = magnitudes.get("SDSS.r")
            bpmag_nc = magnitudes.get("GAIA3.Gbp")

            hpflux_nc = fluxes.get("Hipparcos.Hp_bes")
            btflux_nc = fluxes.get("TYCHO.B_bes")
            vtflux_nc = fluxes.get("TYCHO.V_bes")
            gflux_nc  = fluxes.get("SDSS.g")
            rflux_nc  = fluxes.get("SDSS.r")
            bpflux_nc = fluxes.get("GAIA3.Gbp")


            adp_dati['hpmag_nc ']['value']   = hpmag_nc
            adp_dati['btmag_nc ']['value']   = btmag_nc
            adp_dati['vtmag_nc ']['value']   = vtmag_nc
            adp_dati['gmag_nc ']['value']    = gmag_nc
            adp_dati['rmag_nc ']['value']    = rmag_nc
            adp_dati['bpmag_nc ']['value']   = bpmag_nc
            adp_dati['hpflux_nc ']['value']  = hpflux_nc
            adp_dati['btflux_nc ']['value']  = btflux_nc
            adp_dati['vtflux_nc ']['value']  = vtflux_nc
            adp_dati['gflux_nc ']['value']   = gflux_nc
            adp_dati['rflux_nc ']['value']   = rflux_nc 
            adp_dati['bpflux_nc ']['value']  = bpflux_nc

            spectrum_flat, rms_flat, skew_flat, kurt_flat = flatten_spectrum(wl_, flux_, 300.0)
            
            #spectra PSD   
            #logfreq, psd, psd_smoothed, psd_slope, psd_intercept, psd_slope_e, psd_intercept_e, psd_chi2, psd_r2 = perform_fft(wl_,spectrum_flat)
            logfreq, psd, psd_smoothed, psd_slope, psd_intercept, psd_slope_e, psd_intercept_e, psd_chi2, psd_r2, fft_peak_wl, fft_peak_value = perform_fft(wl_,flux_)
            adp_dati['psd_slope']['value']       = psd_slope
            adp_dati['psd_intercept']['value']   = psd_intercept
            adp_dati['psd_slope_e']['value']     = psd_slope_e
            adp_dati['psd_intercept_e']['value'] = psd_intercept_e
            adp_dati['psd_chi2']['value']        = psd_chi2
            adp_dati['psd_r2']['value']          = psd_r2
            adp_dati['fft_peak_wl']['value']     = fft_peak_wl
            adp_dati['fft_peak_value']['value']  = fft_peak_value


            #adp_dati['']['value']  = 
            #adp_dati['']['value']  = 
            #adp_dati['']['value']  = 


            #update dictionary with results
            all_results.update(adp_data)
            all_results.update(adp_dati)
            

            stat_spec_b = {'max_b': None, 'min_b': None, 'sum_b': None, 'mean_b': None, 'median_b': None, 'rms_b': None, 'skew_b': None, 'kurt_b': None, 'max_bf': None, 'min_bf': None, 'sum_bf': None, 'mean_bf': None, 'median_bf': None, 'rms_bf': None, 'skew_bf': None, 'kurt_bf': None}
            stat_spec_r = {'max_r': None, 'min_r': None, 'sum_r': None, 'mean_r': None, 'median_r': None, 'rms_r': None, 'skew_r': None, 'kurt_r': None, 'max_rf': None, 'min_rf': None, 'sum_rf': None, 'mean_rf': None, 'median_rf': None, 'rms_rf': None, 'skew_rf': None, 'kurt_rf': None}
            
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
                
                if j1 > 0 and j1 < j2:
                    max_b, min_b, sum_b, mean_b, median_b, rms_b, skew_b, kurt_b, rms_bf, min_bf, sum_bf, mean_bf, median_bf, rms_bf, skew_bf, kurt_bf = process_spectrum(flux, wl, j1, j2)
                    stat_spec_b = {'max_b': max_b, 'min_b': min_b, 'sum_b': sum_b, 'mean_b': mean_b, 'median_b': median_b, 'rms_b': rms_b,'skew_b': skew_b, 'kurt_b': kurt_b, 
                                   'max_bf': rms_bf, 'min_bf': min_bf, 'sum_bf': sum_bf, 'mean_bf': mean_bf, 'median_bf': median_bf, 'rms_bf': rms_bf, 'skew_bf': skew_bf, 'kurt_bf': kurt_bf}

                if j3 > j2 and j3 < j4:
                    max_r, min_r, sum_r, mean_r, median_r, rms_r, skew_r, kurt_r, rms_rf, min_rf, sum_rf, mean_rf, median_rf, rms_rf, skew_rf, kurt_rf = process_spectrum(flux, wl, j3, j4)
                    stat_spec_r = {'max_r': max_r, 'min_r': min_r, 'sum_r': sum_r, 'mean_r': mean_r, 'median_r': median_r, 'rms_r': rms_r,'skew_r': skew_r, 'kurt_r': kurt_r, 
                                   'max_rf': rms_rf, 'min_rf': min_rf, 'sum_rf': sum_rf,'mean_rf': mean_rf, 'median_rf': median_rf, 'rms_rf': rms_rf, 'skew_rf': skew_rf, 'kurt_rf': kurt_rf}
            
            all_results.update(stat_spec_b)
            all_results.update(stat_spec_r)


            # RV from H-alpha line
            ha_bw = 20.0
            prefix = 'halpha_'
            values_fit_line = fit_absorption_line(wl_, flux_, ha_wl, ha_bw, True)
            halpha_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
            all_results.update(halpha_results)
            #print(ha_results)

            # correct wavelenght for radial velocity
            wl_rv_corr = wl_ * (1 - halpha_results['halpha_rv'] / ligth_speed)



            fitta_linee = True
            if fitta_linee == True:
                for idx, row in nist_selected_lines.iterrows():
                    prefix = row['Element'] + convert_arabic_to_uppercase_roman_numbers(row['sp_num']) + "_"
                    bw_line = 5
                    wl_line_center = row['lambda_air']
                    values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line, False)
                    linename_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                    all_results.update(linename_results)
        






            other_lines = False
            if other_lines == True:
            # wl_important_lines
                # (4861.332, 10, 'Hbeta'),
                # (4340.468, 10, 'Hgamma'),
                # (4101.737, 10, 'Hdelta'),
                # (3970.074, 10, 'Hepsilon'),
                # (4471.477, 10, 'HeI'),
                # (4541.590, 10, 'HeII'),
    
                #analysis of other prominent lines
                #prefix = 'LINEA_'
                #bw_line = 5
                #wl_line_center = 
                #values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                #LINEA_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                #all_results.update(LINEA_results)
                #print(LINEA_results)
    
                prefix = 'Hbeta_'
                bw_line = 5
                wl_line_center = 4861.332
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                hbeta_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(hbeta_results)
    
                prefix = 'Hgamma_'
                bw_line = 5
                wl_line_center = 4340.468
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                hgamma_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(hgamma_results)
                #print(hgamma_results)
    
                prefix = 'Hdelta_'
                bw_line = 5
                wl_line_center = 4101.737
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                hdelta_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(hdelta_results)
                #print(hdelta_results)
    
                prefix = 'Hepsilon_'
                bw_line = 5
                wl_line_center = 3970.074
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                hepsilon_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(hepsilon_results)
                #print(hepsilon_results)
    
                prefix = 'HeI_'
                bw_line = 5
                wl_line_center = 4471.477
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                hei_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(hei_results)
                #print(hei_results)
    
                prefix = 'HeII_'
                bw_line = 5
                wl_line_center = 4541.590
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                heii_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(heii_results)
                #print(heii_results)
    
    
        # (3933.664, 10, 'CaH'),
        # (4481.129, 10, 'MgII'),
        # (5269.541, 10, 'FeI'),
        # (5895.923, 10, 'NaI_D1'),
        # (5889.953, 10, 'NaI_D2'),
    
                prefix = 'CaH_'
                bw_line = 5
                wl_line_center = 3933.664
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                cah_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(cah_results)
                #print(cah_results)
    
                prefix = 'MgII_'
                bw_line = 5
                wl_line_center = 4481.129
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                mgii_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(mgii_results)
                #print(mgii_results)
    
                prefix = 'FeI_'
                bw_line = 5
                wl_line_center = 5269.541
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                fei_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(fei_results)
                #print(fei_results)
    
                prefix = 'NaI_D1_'
                bw_line = 5
                wl_line_center = 5895.923
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                naid1_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(naid1_results)
                #print(naid1_results)
    
                prefix = 'NaI_D2_'
                bw_line = 5
                wl_line_center = 5889.953
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                naid2_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(naid2_results)
                #print(naid2_results)
    
    
        # (6707.815, 10, 'LiI'),
        # (4554.033, 10, 'BaII'),
        # (4262.270, 10, 'TcI')
    
                prefix = 'LiI_'
                bw_line = 5
                wl_line_center = 6707.815
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                lii_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(lii_results)
                #print(lii_results)
    
    
                prefix = 'BaII_'
                bw_line = 5
                wl_line_center = 4554.033
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                baii_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(baii_results)
                #print(baii_results)
    
    
                prefix = 'TcI_'
                bw_line = 5
                wl_line_center = 4262.270
                values_fit_line = fit_absorption_line(wl_rv_corr, flux_, wl_line_center, bw_line)
                tci_results = {f"{prefix}{var}": value for var, value in zip(variables_fit_line, values_fit_line)}
                all_results.update(tci_results)
                #print(tci_results)






            # spectral classification
            # extract from the spectra only the interval around h-alpha line
            mask_wl_sp_class = (wl_rv_corr > (ha_wl - ha_bw/2)) & (wl_rv_corr < (ha_wl + ha_bw/2))
            wl_spt_class_0   = wl_rv_corr[mask_wl_sp_class]
            flux_sp_class_0  = flux_[mask_wl_sp_class]
            max_flux_sp_class_0 = np.max(flux_sp_class_0)
            flux_sp_class_0 = flux_sp_class_0 / max_flux_sp_class_0

            # extract the wavelenght of the template around h-alpha line
            mask_wl_pha_class = (wl_pha > (ha_wl - ha_bw/2)) & (wl_pha < (ha_wl + ha_bw/2))
            wl_pha_class      = wl_pha[mask_wl_pha_class]
            
            # Further mask wl_pha_class to be within the range of wl_spt_class_0
            min_wl = min(wl_spt_class_0)
            max_wl = max(wl_spt_class_0)
            mask_within_range = (wl_pha_class >= min_wl) & (wl_pha_class <= max_wl)
            wl_pha_class = wl_pha_class[mask_within_range]


            # interpolate the spectra on the template wavelegth
            fl_interpolator = interp1d(wl_spt_class_0, flux_sp_class_0, kind='linear')
            flux_sp_class = fl_interpolator(wl_pha_class)
            
            # find the best matching template
            chi2_class = 9e9
            i_class = -1
            classification = defaultdict(dict)
            for ii in range(len(pha_selected)):
                spec_pha_0 = np.array(pha_selected['spectrum_slice'][ii])
                spec_pha = spec_pha_0[mask_wl_pha_class]
                spec_pha = spec_pha[mask_within_range]/np.max(spec_pha[mask_within_range])
                chi2 = np.sum((flux_sp_class - spec_pha)**2)
                if chi2 < chi2_class:
                    chi2_class = chi2
                    i_class = ii
            if i_class>0:
                classification['mass_class']['value']   = pha_selected['Mass'][i_class]
                classification['radius_class']['value'] = pha_selected['Radius'][i_class]
                classification['logg_class']['value']   = pha_selected['logg'][i_class]
                classification['meh_class']['value']    = pha_selected['MeH'][i_class]
                classification['logL_class']['value']   = pha_selected['LogL'][i_class]
                classification['teff_class']['value']   = 10**pha_selected['logTe'][i_class]
                classification['logrho_class']['value'] = pha_selected['logRho'][i_class]
                classification['chi2_class']['value'] = chi2
            else:
                classification['mass_class']['value']   = np.nan
                classification['radius_class']['value'] = np.nan
                classification['logg_class']['value']   = np.nan
                classification['meh_class']['value']    = np.nan
                classification['logl_class']['value']   = np.nan
                classification['teff_class']['value']   = np.nan
                classification['logrho_class']['value'] = np.nan
                classification['chi2_class']['value']   = np.nan
            #print(classification)

            # analyze the spectral curvature
            wl_rebin,spectrum_rebin = rebin_spectrum(wl_, flux_, 500.0)
            result_curvature = analyze_curvature(wl_rebin, spectrum_rebin)
            classification['curvature']['value'] = result_curvature
            all_results.update(classification)
            
            

            harps_obs.close()

# ++++++++++++++++++++++++  END   ADP SPECTRA  +++++++++++++++++++++++ #


        # Create patterns to search for "bis", "ccf", and "guide" files
        # Check if the length of archive_id_ancillary is greater than 0 and if there is at least one "bis" or one "ccf" or one "guide" file

# ++++++++++++++++++++++++  START CCF FILE     +++++++++++++++++++++++ #
        ccf_pattern = os.path.join(storage_path, dateobs, dp_id_raw + "_ccf*_A.fits")
        ccf_files = glob.glob(ccf_pattern)
        ccf_data = process_file(ccf_files, 'ccf', datakw)
        all_results.update(ccf_data)
# ++++++++++++++++++++++++  END CCF FILE       +++++++++++++++++++++++ #

# ++++++++++++++++++++++++  START BIS FILE     +++++++++++++++++++++++ #
        bis_pattern = os.path.join(storage_path, dateobs, dp_id_raw + "_bis*_A.fits")
        bis_files = glob.glob(bis_pattern)
        bis_data = process_file(bis_files, 'bis', datakw)
        all_results.update(bis_data)
# ++++++++++++++++++++++++  END BISECTOR FILE  +++++++++++++++++++++++ #

# ++++++++++++++++++++++++  START GUIDE FILE   +++++++++++++++++++++++ #
        gui_pattern = os.path.join(storage_path, dateobs, dp_id_raw + "_INT_GUIDE.fits")
        gui_files = glob.glob(gui_pattern)
        gui_data = process_gui_file(gui_files, datakw)
        all_results.update(gui_data)
# ++++++++++++++++++++++++  END GUIDE FILE     +++++++++++++++++++++++ #
        

        all_results.update(filenames_dict)

        # print all the keys in the dictionary
        #print(len(all_results))
        #keys_list = list(all_results.keys())  # Convert dictionary keys to a list
        #keys_str = ', '.join(keys_list)  # Join the keys with commas
        #print(keys_str)

        out_values = []
        for value in all_results.values():
            if isinstance(value, dict):
                out_values.append(value['value'])
            else:
                out_values.append(value)


        if out_values:
            with open(outfile, 'a',newline='') as csv_file:
                writer = csv.writer(csv_file, lineterminator='\n')
                writer.writerow(out_values)








sys.exit()
#---------------   END PROGRAM









            

            #trasm_obs=interpolate_spectrum(harps_wavelength0,harps_trasm_all,wl_)
            #fig, ax1 = plt.subplots(1, 1)
            #ax1.plot(wl_,flux_/np.mean(flux_))
            #ax1.plot(wl_,trasm_obs*10)
            #ax1.plot(wl_,flux_/trasm_obs)
            #plt.show()

            
            #print(adp_data['object']['value'])
            # Tecnezio
            #popt,deltawl,rvha = fit_absorption_line(wl_, flux_,4238.19,3.0)
            #popt,deltawl,rvha = fit_absorption_line(wl_, flux_,4262.27,3.0)
            #popt,deltawl,rvha = fit_absorption_line(wl_, flux_,4297.00,3.0)
            
            #print()
            #print(row['object'],row['spt_obj'])
            #rv,reduced_chisq,gamma,skewness,kurtosis = fit_absorption_line2(wl_, flux_,6562.817,30.0)
            
            
            #if k==1:
            #    speck = np.zeros_like(speck1)

            #print(k)
            #speck=((k-1)*speck+speck1)/k
            
            #if np.mod(k,100)==0:
            #    fig, ax1 = plt.subplots(1, 1)
            #    ax1.plot(speck)
            #    plt.show()
            

            
            #new_wl = wl_*(1+rv/const.c.to(u.km/u.s).value)
            #flux_interpolator = interp1d(wl_, flux_, kind='cubic', bounds_error=False, fill_value=0)
            #new_flux = flux_interpolator(new_wl)
            #_,_,_,_,_ = fit_absorption_line(new_wl, new_flux,6562.817,10.0, True)

            
            #popt,deltawl,rvha = fit_absorption_line(wl_, flux_,4861.332,10.0)
            #continuum,intensity,centrawl,width1sig = popt
            #print(row['object'],",",row['mjd_obs'],",",rv,",",reduced_chisq,",",gamma,",",skewness,",",kurtosis)








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

# ------------------------------  BEGIN FUNCTION ------------------------------ #
# def fit_absorption_line(wavelengths, intensities, wc,ws, flagd):
    # plot=flagd
    # log=flagd
    # line_center0 = wc
    # wl_span = ws
    # mask = (wavelengths > (line_center0 - wl_span)) & (wavelengths < (line_center0 + wl_span))
    # wl=wavelengths[mask]
    # flux = intensities[mask]
    # npt=len(wavelengths[mask])
    
    # #for attempt in range(5):
    # const_guess = np.mean(flux)
    # amplitude_guess = np.min(flux) - const_guess
    # mean_guess = line_center0 #*(attempt+1)
    # stddev_guess = 0.1

    # #fit Gaussian
    # ERRORE=0
    # try:
        # popt_gaussian, pcov = curve_fit(negative_gaussian, wl, flux, p0=[const_guess, amplitude_guess, mean_guess, stddev_guess], maxfev=2000)
    # except RuntimeError as e:
        # ERORRE = 1
        # popt_gaussian= [np.nan for _ in range(4)]
        # pcov = [[np.nan for _ in range(4)] for _ in range(4)]
        # errors = [np.nan for _ in range(4)]

    # if ERRORE == 0:
        # y_fit_gaussian = negative_gaussian(wl, *popt_gaussian)
        # errors = np.sqrt(np.diag(pcov))
        # const_c=popt_gaussian[0]
        # depth_c=popt_gaussian[1]
        # lambda_c=popt_gaussian[2]
        # sigma=popt_gaussian[3]
        # sigma_error = errors[3]

    # # Fit Lorentzian
    # ERRORE=0
    # try:
        # if np.abs(const_c)>0 and np.abs(depth_c)>0 and np.abs(lambda_c)>0:
            # const_guess=const_c
            # amplitude_guess=depth_c
            # mean_guess=lambda_c
            # stddev_guess=sigma
        # popt_lorentzian, pcov = curve_fit(negative_lorentzian, wl, flux, p0=[const_guess, amplitude_guess, mean_guess, stddev_guess], maxfev=2000)
    # except RuntimeError as e:
        # ERORRE = 1
        # popt_lorentzian= [np.nan for _ in range(4)]
        # pcov = [[np.nan for _ in range(4)] for _ in range(4)]
        # errors = [np.nan for _ in range(4)]
        # delta_wl=np.nan
        # rv=np.nan

    # if ERRORE == 0:
        # y_fit_lorentzian = negative_lorentzian(wl, *popt_lorentzian)
        # errors = np.sqrt(np.diag(pcov))  # Extract the errors from the covariance matrix
    
        # continuum=popt_lorentzian[0]
        # amplitude=popt_lorentzian[1]
        # line_center=popt_lorentzian[2]
        # gamma=popt_lorentzian[3]
        # continuum_error = errors[0]
        # amplitude_error = errors[1]
        # line_center_error = errors[2]
        # gamma_error = errors[3]
        # delta_wl=line_center0-popt_lorentzian[2]
        # rv = (line_center/line_center0 - 1)*const.c.to(u.km/u.s).value
        
        # # Subset of flux and y_fit
        # mask1 = (wl > (line_center - 2 * gamma)) & (wl < (line_center + 2 * gamma))
        # wl_subset=wl[mask1]
        # flux_subset = flux[mask1]/continuum
        # y_fit_lorentzian_subset = y_fit_lorentzian[mask1]/continuum
        
        # residuals = flux_subset - y_fit_lorentzian_subset
        # rms_residuals = np.std(residuals)
        # snr=abs(amplitude/continuum)/rms_residuals
        
        # chisq = np.sum((residuals / np.std(flux_subset)) ** 2)
        # dof = len(wl) - len(popt_lorentzian)
        # reduced_chisq = chisq / dof
        # ss_res = np.sum(residuals ** 2)
        # ss_tot = np.sum((flux_subset - np.mean(flux_subset)) ** 2)
        # r_squared = 1 - (ss_res / ss_tot)
        # if r_squared<0 :
            # print("************************************")
            # print("************************************")
            # print("* W A R N I N G       R2<0          ")
            # print("************************************")
            # print("************************************")
            # #    break
            
    
    # #if ERRORE == 0:
        # #print('attempt ',attempt)
        # #if attempt==1000:
        # #        plot=False
        # #        log=False
                
        # flux1=1-flux_subset
        # sum_flux = sum(flux1)
        # npt_subset = len(flux_subset)
        
        # skewness = sum(flux1[i] * (wl_subset[i] - line_center)**3 for i in range(npt_subset)) / (sum_flux * sigma**3)
        # kurtosis = sum(flux1[i] * (wl_subset[i] - line_center)**4 for i in range(npt_subset)) / (sum_flux * sigma**4) - 3

        # #if amplitude>0 and line_center0>line_center-gamma and line_center0<line_center+gamma and abs(line_center-line_center0)<0.25:
        # if plot == True:    
    # #       log=True
            # if log == True:
                # print('continuum         ',continuum, continuum_error)
                # print('intensity         ',amplitude, amplitude_error)
                # print('central wavelength',line_center,line_center_error)
                # print('gamma   = width   ',gamma, gamma_error)
                # print('sigma   = width   ',sigma)
                # print('skewness          ',skewness)
                # print('kurtosis          ',kurtosis)
                # print('rms residuals     ',rms_residuals)
                # print('SNR               ',snr)
                # print('chi2              ',reduced_chisq)
                # print('R2                ',r_squared)
                # print('delta wavelength  ',delta_wl)
                # print('RV                ',rv)
            # #Blue: '#1F77B4'
            # #Orange: '#FF7F0E'
            # #Green: '#2CA02C'
            # #Red: '#D62728'
            # #Purple: '#9467BD'
    # #        plot = True
     
            # if plot == True:
                # print(rv/const.c.to(u.km/u.s).value)
                # x1 = line_center0 - wl_span
                # x2 = line_center0 + wl_span
                # #x1 = line_center0 - 2
                # #x2 = line_center0 + 2
                # lambda_values = absorption_lines_data['lambda_air']
                # loggf_values = absorption_lines_data['log_gf']
                # elements = absorption_lines_data['specie']
                # mask_lines = (lambda_values > (line_center0 - wl_span)) & (lambda_values < (line_center0 + wl_span)) 
                # #mask_lines2 = (lambda_values > (line_center0 - wl_span)) & (lambda_values < (line_center0 + wl_span)) & np.isnan(loggf_values)
                # #new_wl = wl*(1+rv/const.c.to(u.km/u.s).value)
                # #flux_interpolator = interp1d(new_wl, flux, kind='cubic', bounds_error=False, fill_value=0)
                # #new_flux = flux_interpolator(wl)
                # new_wl = wl*(1-rv/const.c.to(u.km/u.s).value)
                # new_flux = flux
                # print("Radial velocity (km/s):", rv)
                # print("Speed of light (km/s):", const.c.to(u.km/u.s).value)
                # print("Expected shift in Angstroms:", 6562 * (rv / const.c.to(u.km/u.s).value))
                # print('mean shift in wavelength',np.mean(wl-new_wl),'mean difference in flux',np.mean(flux-new_flux))
                # #for w, f, new_w, new_f in zip(wl, flux, new_wl, new_flux):
                # #    print(w, f, new_w, new_f)
                # print()
                # num_plot = 1
                # if num_plot == 1:
                    # fig, ax1 = plt.subplots(1, 1)
                    # manager = plt.get_current_fig_manager()
                    # manager.window.wm_geometry("1400x1200+0+0")
                    # original = False
                    # if original == True:
                        # ax1.plot(new_wl,new_flux/continuum, color='orange', linestyle ='dashed', label='Data')
                        # ax1.axvline(line_center0, color='#00FF00', linestyle='dashed', label='Center')
                        # ax1.axvline(line_center , color='#FF00FF', linestyle='dashed', label='Center')
                        # ax1.axvline(line_center-gamma, color='#FF00FF', linestyle='dashed', label='Center')
                        # ax1.axvline(line_center+gamma, color='#FF00CC', linestyle='dashed', label='Center')
                    # else:
                        # ax1.plot(new_wl,new_flux/continuum, label='Data')
                        # ax1.plot(new_wl, y_fit_lorentzian/continuum, label='Lorentzian Fit', linestyle='dashed')
                        # ax1.axvline(line_center0, color='#00FF00', linestyle='dashed', label='Center')
                        # jj=-1
                        # for lam, elem in zip(lambda_values[mask_lines], elements[mask_lines]):
                            # jj= jj+1
                            # ax1.axvline(lam, color = '#bbbbff', linestyle='dotted')
                            # ax1.text(lam, ax1.get_ylim()[0]+0.0001*jj, elem, rotation=45, va='bottom')
                        # #for lam, elem in zip(lambda_values[mask_lines2], elements[mask_lines2]):
                        # #    ax1.axvline(lam, linestyle='dotted')
                        # #    ax1.text(lam, ax1.get_ylim()[0], elem, rotation=45, va='bottom')
                    # ax1.set_xlim(x1, x2)
                    # ax1.set_ylabel('Flux')
                    # ax1.grid(True)
                    # plt.tight_layout()
                    # plt.show()
                # elif num_plot == 2:
                    # fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
                    # manager = plt.get_current_fig_manager()
                    # manager.window.wm_geometry("1800x600+0+0")
                    # original = False
                    # if original == True:
                        # ax1.plot(new_wl,new_flux/continuum, color='orange', linestyle ='dashed', label='Data')
                        # ax1.axvline(line_center0, color='#00FF00', linestyle='dashed', label='Center')
                        # ax1.axvline(line_center , color='#FF00FF', linestyle='dashed', label='Center')
                        # ax1.axvline(line_center-gamma, color='#FF00FF', linestyle='dashed', label='Center')
                        # ax1.axvline(line_center+gamma, color='#FF00CC', linestyle='dashed', label='Center')
                    # else:
    # #                    ax1.plot(wl,flux/continuum, color ='#fdfdfd', label='Data')
                        # ax1.plot(new_wl,new_flux/continuum, label='Data')
                        # ax1.plot(new_wl, y_fit_lorentzian/continuum, label='Lorentzian Fit', linestyle='dashed')
                        # ax1.axvline(line_center0, color='#00FF00', linestyle='dashed', label='Center')
                        # for lam, elem in zip(lambda_values[mask_lines], elements[mask_lines]):
                            # ax1.axvline(lam, color='#00e000', linestyle='dotted')
                            # ax1.text(lam, ax1.get_ylim()[0], elem, rotation=45, va='top', color='#00a000')
                        # for lam, elem in zip(lambda_values[mask_lines2], elements[mask_lines2]):
                            # ax1.axvline(lam, color='#d000b0', linestyle='dotted')
                            # ax1.text(lam, ax1.get_ylim()[0], elem, rotation=45, va='bottom', color='#8000b0')
                    # ax1.set_xlim(x1, x2)
                    # ax1.set_ylabel('Flux')
                    # #ax1.set_title('H-alpha Line')
                    # ax1.grid(True)
                    # #ax1.legend()
                    # ax2.plot(new_wl, (flux - y_fit_lorentzian)/continuum, label='Residuals')
                    # ax2.axvline(line_center0, color='#00FF00', linestyle='dashed', label='Center')
                    # #ax2.axvline(line_center-gamma, color='red', linestyle='dashed', label='Center')
                    # #ax2.axvline(line_center+gamma, color='red', linestyle='dashed', label='Center')
                    # for lam, elem in zip(lambda_values[mask_lines], elements[mask_lines]):
                        # ax2.axvline(lam, color='#000080', linestyle='dotted')
                    # ax2.set_xlim(x1, x2)
                    # ax2.set_xlabel('Wavelength [AA]')
                    # ax2.set_ylabel('Residuals')
                    # ax2.grid(True)
                    # plt.tight_layout()
                    # plt.show()
                # else:
                    # figatomare = True
            
        # new_wl = wl*(1-rv/const.c.to(u.km/u.s).value)
        # new_flux = flux
        # speck=new_flux/continuum
        # speck_wl=new_wl
    # else:
        # figatomare = True

    # return rv,reduced_chisq,gamma,skewness,kurtosis,speck,speck_wl
# # ------------------------------  END   FUNCTION ------------------------------ #
