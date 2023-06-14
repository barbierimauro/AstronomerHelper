#!/usr/bin/python3
import sys
import os
import pandas as pd
import numpy as np
import scipy.stats as st
from astropy.io import fits
from astropy.table import Table
from scipy.ndimage import uniform_filter1d
from scipy.interpolate import interp1d
from scipy.integrate import simpson
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def planck(wavelength, temperature):
    numerator = 2 * plank_constant * ligth_speed ** 2
    exponent = (plank_constant * ligth_speed) / (wavelength * boltzmann_constant * temperature)
    denominator = (wavelength ** 5) * (np.exp(exponent) - 1)
    black_body = numerator / denominator
    return black_body

def safe_value(value, placeholder=-9999):
    if np.isnan(value) or np.isinf(value):
        return placeholder
    else:
        return value

def mag_integral(y_values, x_values):
    if len(x_values) % 2 == 0:
        x_values = x_values[:-1]
        y_values = y_values[:-1]
    integral = simpson(y_values*x_values, x_values)
    #integral = simpson(y_values, x_values)
    #integral = np.sum(y_values)
    magnitude = -2.5*np.log10(integral)
    return magnitude


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


##########################################################################
# START MAIN CODE

# CONSTANTS
solar_mass=1.9885e33 # g
solar_radius=6.957e10 #cm
solar_luminosity=3.828e33 # erc/cm**2/s
plank_constant = 6.62607004e-34  # Planck's constant (J s)
ligth_speed = 3.0e8  # Speed of light (m/s)
boltzmann_constant = 1.38064852e-23  # Boltzmann's constant (J/K)
wielen_constant = 2.897771955e7
au = 14959787070000.0 # astronomical unit in cm
pc = 3085677581491367300.0 # parsec in cm

# DEFINITIONS
window_size = 100          #number of points
wlmin=3700                 # Angstrom
wlmax=7000                 # Angstrom
#wlmin=3000                 # Angstrom
#wlmax=53000                # Angstrom
wl_min_b=3783              # Angstrom
wl_max_b=5304              # Angstrom
wl_min_r=5338              # Angstrom
wl_max_r=6912              # Angstrom
coarse_resolution = 100.0  # Angstrom

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
harps_trasm_all0=harps_trasm_det0*harps_trasm_fib0*harps_trasm_inst0*harps_trasm_tel0

#wavelength_range       = (harps_wavelength0 >= wl_min_b) & (harps_wavelength0 <= wl_max_b)

selected_wl_b_harps    = harps_wavelength0[(harps_wavelength0 >= wl_min_b) & (harps_wavelength0 <= wl_max_b)]
selected_trasm_b_harps = harps_trasm_all0[(harps_wavelength0  >= wl_min_b) & (harps_wavelength0 <= wl_max_b)]
selected_atm_b_harps   = harps_trasm_atm0[(harps_wavelength0  >= wl_min_b) & (harps_wavelength0 <= wl_max_b)]
selected_bla_b_harps   = harps_trasm_bla0[(harps_wavelength0  >= wl_min_b) & (harps_wavelength0 <= wl_max_b)]
selected_dis_b_harps   = harps_dispersion0[(harps_wavelength0 >= wl_min_b) & (harps_wavelength0 <= wl_max_b)]

selected_wl_r_harps    = harps_wavelength0[(harps_wavelength0 >= wl_min_r) & (harps_wavelength0 <= wl_max_r)]
selected_trasm_r_harps = harps_trasm_all0[(harps_wavelength0  >= wl_min_r) & (harps_wavelength0 <= wl_max_r)]
selected_atm_r_harps   = harps_trasm_atm0[(harps_wavelength0  >= wl_min_r) & (harps_wavelength0 <= wl_max_r)]
selected_bla_r_harps   = harps_trasm_bla0[(harps_wavelength0  >= wl_min_r) & (harps_wavelength0 <= wl_max_r)]
selected_dis_r_harps   = harps_dispersion0[(harps_wavelength0 >= wl_min_r) & (harps_wavelength0 <= wl_max_r)]







csv_file = 'phoenix_spectral_library_hires.csv'
df = pd.read_csv(csv_file)

hdu_wl = fits.open('WAVE_PHOENIX-ACES-AGSS-COND-2011.fits', ignore_missing_end=True)
wl = hdu_wl[0].data
hdu_wl.close()

wavelength_range = (wl >= wlmin) & (wl <= wlmax)
selected_wl = wl[wavelength_range]
w_b = selected_wl[(selected_wl >= wl_min_b) & (selected_wl <= wl_max_b)]
w_r = selected_wl[(selected_wl >= wl_min_r) & (selected_wl <= wl_max_r)]
w_h = selected_wl



#sys.exit()

output_csv = 'phoenix_harps_spectra.csv'
header="FILENAME,MASS,RADIUS,LOGG,MEH,LOGL,LOGTE,LOGRHO,PEAKBB,BOLMAG,HMAG,BMAG,RMAG,BRCOL,HBCOL,HRCOL,MAX_B,MIN_B,SUM_B,MEAN_B,MEDIAN_B,RMS_B,SKEW_B,KURT_B,MAX_R,MIN_R,SUM_R,MEAN_R,MEDIAN_R,RMS_R,SKEW_R,KURT_R,MAX_H,MIN_H,SUM_H,MEAN_H,MEDIAN_H,RMS_H,SKEW_H,KURT_H,RMSFLAT,SKEWFLAT,KURTFLAT"
if os.path.exists(output_csv):
    os.remove(output_csv)

if not os.path.exists(output_csv):
    with open(output_csv, 'w') as f:
        f.write(header + "\n")

vega_file='PHOENIX-ACES-AGSS-COND-2011/Z-0.5/lte09600-4.00-0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
vega_hdu = fits.open(vega_file, ignore_missing_end=True)
vega_hdr = vega_hdu[0].header
vega_spec_full = vega_hdu[0].data
vega_hdu.close()

vega_spec_selected = vega_spec_full[wavelength_range]
vega_radius = vega_hdr["PHXREFF"]
vega_rconst = -2.5*np.log10((vega_radius/ (10*pc))**2)

vega_mbol = mag_integral(vega_spec_full, wl) + vega_rconst
x_0 = vega_spec_selected[(selected_wl >= wl_min_b) & (selected_wl <= wl_max_b)]
vega_bmag = mag_integral(x_0, w_b) + vega_rconst
x_0 = vega_spec_selected[(selected_wl >= wl_min_r) & (selected_wl <= wl_max_r)]
vega_rmag = mag_integral(x_0, w_r) + vega_rconst
x_0 = vega_spec_selected
vega_hmag = mag_integral(x_0, w_h) + vega_rconst


print('vega',vega_bmag,vega_rmag,vega_hmag,vega_mbol)
print()
vega_mbol = mag_integral(vega_spec_full, wl) + vega_rconst
x_0 = vega_spec_selected[(selected_wl >= wl_min_b) & (selected_wl <= wl_max_b)]
vega_bmag = mag_integral(x_0 * selected_trasm_b_harps, w_b) + vega_rconst
x_0 = vega_spec_selected[(selected_wl >= wl_min_r) & (selected_wl <= wl_max_r)]
vega_rmag = mag_integral(x_0 * selected_trasm_r_harps, w_r) + vega_rconst
x_0 = vega_spec_selected
vega_hmag = mag_integral(x_0, w_h) + vega_rconst


print('vega',vega_bmag,vega_rmag,vega_hmag,vega_mbol)
print()
sys.exit()


for index, row in df.iterrows():
    file_name = row['fname']
    file_path = 'PHOENIX-ACES-AGSS-COND-2011/' + file_name
    # test with an A0I star
    #file_path = 'PHOENIX-ACES-AGSS-COND-2011/'+'Z+0.5/lte10000-2.00+0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
    hdu = fits.open(file_path, ignore_missing_end=True)
    hdr = hdu[0].header
    spec = hdu[0].data
    abundances_h = hdu[1].header
    abundances = hdu[1].data
    hdu.close()

    ##############################################
    #

    selected_spec = spec[wavelength_range]
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #####SCOMMENTARE
    #selected_err = rolling_std(selected_spec, window_size)
    selected_err = selected_spec
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    radius = hdr["PHXREFF"]
    rconst = -2.5*np.log10((radius/ (10*pc))**2)

    bolmag = mag_integral(spec, wl) - vega_mbol + rconst
    
    x_0 = selected_spec[(selected_wl >= wl_min_b) & (selected_wl <= wl_max_b)]
    max_b, min_b, sum_b, mean_b, median_b, rms_b, skew_b, kurt_b = statistics(x_0)
    bmag = mag_integral(x_0, w_b) - vega_bmag + rconst

    x_0 = selected_spec[(selected_wl >= wl_min_r) & (selected_wl <= wl_max_r)]
    max_r, min_r, sum_r, mean_r, median_r, rms_r, skew_r, kurt_r = statistics(x_0)
    rmag = mag_integral(x_0, w_r) - vega_rmag + rconst 

    x_0 = selected_spec
    max_h, min_h, sum_h, mean_h, median_h, rms_h, skew_h, kurt_h = statistics(x_0)
    hmag = mag_integral(x_0, w_h) - vega_hmag + rconst


    
    brcol = bmag - rmag
    hbcol = hmag - bmag
    hrcol = hmag - rmag

    print(bmag,rmag,hmag,bolmag)
    print(brcol,hbcol,hrcol)

    
    MASS=hdr["PHXMASS"]/solar_mass
    RADIUS=hdr["PHXREFF"]/solar_radius
    LOGG=hdr["PHXLOGG"]
    MEH=hdr["PHXM_H"]
    LOGL=np.log10(hdr["PHXLUM"]/solar_luminosity)
    LOGTE=np.log10(hdr["PHXTEFF"])
    LOGRHO=np.log10(3/(4*np.pi)*MASS/RADIUS**3)
    
    spectrum_flat, rms_flat, skew_flat, kurt_flat = flatten_spectrum(selected_wl, selected_spec, coarse_resolution)

    wlbb=np.linspace(wlmin*1e-10, wlmax*1e-10, len(selected_wl))
    black_body = planck(wlbb, hdr["PHXTEFF"])
    peakbb = wielen_constant/hdr["PHXTEFF"]
    norm_spec_star=(selected_spec/max(selected_spec))
    norm_spec_black_body=(black_body/max(black_body))

    ##############################################
    #
    
    table_data = Table([selected_wl, selected_spec, selected_err, black_body, spectrum_flat,norm_spec_star,norm_spec_black_body], names=('WAVE', 'FLUX', 'ERR', 'BLACK_BODY','FLAT_SPEC','NORM_SPEC_STAR','NORM_SPEC_BLACK_BODY'))
    new_hdu = fits.BinTableHDU(table_data)

    ##############################################
    #

    primary_hdr = fits.Header()
    for keyword in hdr:
        if keyword.startswith("PH") or keyword in ["BUNIT", "DATE", "GRIDVER", "WAVE"]:
            primary_hdr[keyword] = hdr[keyword]

    primary_hdr.update({
        'MASS': safe_value(MASS),
        'RADIUS': safe_value(RADIUS),
        'LOGG': safe_value(LOGG),
        'MEH': safe_value(MEH),
        'LOGL': safe_value(LOGL),
        'LOGTE': safe_value(LOGTE),
        'LOGRHO': safe_value(LOGRHO),
        'PEAKBB': safe_value(peakbb),
        'BOLMAG': safe_value(bolmag),
        'HMAG': safe_value(hmag),
        'BMAG': safe_value(bmag),
        'RMAG': safe_value(rmag),
        'BRCOL': safe_value(brcol),
        'HBCOL': safe_value(hbcol),
        'HRCOL': safe_value(hrcol),
        'MAX_B': safe_value(max_b),
        'MIN_B': safe_value(min_b),
        'SUM_B': safe_value(sum_b),
        'MEAN_B': safe_value(mean_b),
        'MEDIAN_B': safe_value(median_b),
        'RMS_B': safe_value(rms_b),
        'SKEW_B': safe_value(skew_b),
        'KURT_B': safe_value(kurt_b),
        'MAX_R': safe_value(max_r),
        'MIN_R': safe_value(min_r),
        'SUM_R': safe_value(sum_r),
        'MEAN_R': safe_value(mean_r),
        'MEDIAN_R': safe_value(median_r),
        'RMS_R': safe_value(rms_r),
        'SKEW_R': safe_value(skew_r),
        'KURT_R': safe_value(kurt_r),
        'MAX_H': safe_value(max_h),
        'MIN_H': safe_value(min_h),
        'SUM_H': safe_value(sum_h),
        'MEAN_H': safe_value(mean_h),
        'MEDIAN_H': safe_value(median_h),
        'RMS_H': safe_value(rms_h),
        'SKEW_H': safe_value(skew_h),
        'KURT_H': safe_value(kurt_h),
        'RMSFLAT': safe_value(rms_flat),
        'SKEWFLAT': safe_value(skew_flat),
        'KURTFLAT': safe_value(kurt_flat)
    })

#    keywords = list(primary_hdr.keys())
#    excluded_keywords = ['SIMPLE','BITPIX','NAXIS','EXTEND','BUNIT','WAVE',
#    'DATE','END','PHXTEFF','PHXLOGG','PHXM_H','PHXALPHA','PHXDUST','PHXEOS',
#    'PHXBUILD','PHXVER','PHXXI_L','PHXXI_M','PHXXI_N','PHXMASS','PHXREFF',
#    'PHXLUM','GRIDVER','PHXMXLEN','PHXCONV']
#    keywords = [k for k in keywords if k not in excluded_keywords]
#    keyword_string = ",".join(keywords)
#    print(keyword_string)
#    sys.exit()




    ##############################################
    #

    sec_hdu = fits.BinTableHDU(abundances,abundances_h)
  
    primary_hdu = fits.PrimaryHDU(header=primary_hdr)
    hdul = fits.HDUList([primary_hdu, new_hdu, sec_hdu])
    
    ##############################################
    #

    new_hdu.header.update({
        'EXTNAME': 'SPECTRUM',
        'NELEM': len(selected_wl),
        'TUNIT1': 'Angstrom',
        'TUNIT2': hdr['BUNIT'],
        'TUNIT3': hdr['BUNIT'],
        'TDMIN1': min(selected_wl),
        'TDMAX1': max(selected_wl),
        'TDMIN2': min(selected_spec),
        'TDMAX2': max(selected_spec),
        'TDMIN3': min(selected_err),
        'TDMAX3': max(selected_err),
        'TUTYP1': 'Spectrum.Data.SpectralAxis.Value',
        'TUCD1': 'em.wl;obs.atmos',
        'TUTYP2': 'Spectrum.Data.FluxAxis.Value',
        'TUCD2': 'phot.flux.density;em.wl',
        'TUTYP3': 'Spectrum.Data.FluxAxis.Accuracy.StatError',
        'TUCD3': 'stat.error;phot.flux.density;em.wl'
    })

    ##############################################
    #

    base_name = file_name.split("/")[-1][:-len(".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")]
    outfile = f"phoenix_for_harps/phoenix_harps_{index + 1:04d}_{base_name}.fits"
    #hdul.writeto(outfile, overwrite=True)

    ##############################################
    #

    with open(output_csv, 'a') as f:
        line_data = [file_name]
        for keyword in header.split(",")[1:]:
            if keyword in primary_hdr:
                line_data.append(str(primary_hdr[keyword]))
            else:
                line_data.append('nan')
        line_str = ','.join(line_data)
        f.write(line_str + "\n")    

    print(index, file_name)
    print()
    #sys.exit()



