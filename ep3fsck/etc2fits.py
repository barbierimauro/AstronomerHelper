#!/usr/bin/python3 
# (c) 2023 MAURO BARBIERI   maurobarbieri.science@gmail.com 
""" 

This script is a converter for ESO ETC (European Southern Observatory 
Exposure Time Calculator) output files. The input JSON file is a 
collection of spectral information from the ETC output, which can be 
quite hard to read and work with. This script aims to make the output 
data more accessible by converting it into a human-readable CSV or FITS 
file.

The script first reads the input JSON file and extracts the necessary 
data contained within nested dictionaries and lists. The JSON file 
contains data for different orders, detectors, and wavelengths. The 
script processes this complex structure to obtain the required data 
fields and append them to the output data.

For each data field, the script checks if the corresponding group and 
series exist in the input JSON data. If the data is not available, it 
appends an appropriate default value (None for FITS and an empty string 
for CSV) for the output file format.

The output data is then sorted based on the wavelength and written to 
either a CSV or FITS file, depending on the file extension specified by 
the user. In both cases, the output file contains headers with 
metadata, including the processing software, version number and date.

Usage:
  python script_name.py -i INPUT_JSON -o OUTPUT_FILE

Arguments:
  -i, --input      : Input JSON file containing ESO ETC data
  -o, --output     : Output CSV or FITS file to store the converted data

Example:
  python script_name.py -i input_data.json -o output_data.fits


"""

import json
from pathlib import Path
import argparse
import sys
import os
import time
from astropy.table import Table
from astropy.io import fits
import csv
from datetime import datetime

start_time = time.time()

SWNAME = "ETC2FITS"
AUTHOR = "MAURO BARBIERI"
EMAIL = "maurobarbieri.science@gmail.com"
VERSION = '1.0'
VERSION_DATE = "2023-05-08"

print("ESO ETC 2 FITS/CSV CONVERTER")
print(f"(c) 2023 {AUTHOR} -- {EMAIL}")
print()

parser = argparse.ArgumentParser(description="Write the output of the ESO ETC from NASTY json format to FITS or CSV")
parser.add_argument("-i", "--input", dest="inputfilename", help="Input json file", required=False)
parser.add_argument("-o", "--output", dest="outputfilename", help="Output FITS or CSV file. Allowed extensions: .fits, .csv", required=False)
args, unknown = parser.parse_known_args()

inputfilename = args.inputfilename
outputfilename = args.outputfilename

if unknown:
    parser.print_help()
    print()
    print("the following arguments are required: -i/--input, -o/--output")
    print("the supported output format are only FITS or CSV")
    print()
    sys.exit(1)

if args.inputfilename is None or args.outputfilename is None:
    parser.print_help()
    print()
    print("the following arguments are required: -i/--input, -o/--output")
    print("the supported output format are only FITS or CSV")
    print()
    sys.exit(1)


file_extension = Path(outputfilename).suffix.lower()
if file_extension not in ['.fits', '.csv']:
    print(f"Unsupported file extension '{file_extension}', using FITS as default.")
    file_extension = '.fits'

if not Path(inputfilename).is_file():
    print(f"Input file '{inputfilename}' does not exist.")
    sys.exit(1)

if Path(outputfilename).is_file():
    print(f"Output file '{outputfilename}' already exists.")
    sys.exit(1)

with open(inputfilename) as f:
    data = json.load(f)

headers = [
    "wavelength", "throughput_atmosphere", "throughput_blaze", "throughput_detector", "throughput_fibinj",
    "throughput_instrument", "throughput_telescope", "throughput_totalinclsky", "dispersion_dispersion", "sed_sky",
    "sed_target", "signals_obssky", "signals_obstarget", "signals_obstotal", "snr_snr", "maxsignals_maxpixelsky",
    "maxsignals_maxpixeltarget", "maxsignals_maxpixeltotal", "Order"
]

output_data = []

# Iterate through the orders in the JSON data
for order_data in data['data']['orders']:
    # Get the order number
    order_number = order_data['order']
    # Iterate through the detectors in the order
    for detector_data in order_data['detectors']:
        # Get the detector series data
        detector_series = detector_data['data']
        # Extract and scale wavelengths
        wavelengths = [round(w * 1e10, 2) for w in detector_series['wavelength']['wavelength']['data']]
        # Iterate through the wavelengths
        for i, wavelength in enumerate(wavelengths):
            # Start a new row with the wavelength
            row_data = [wavelength]
            # Iterate through the headers (excluding the first and last headers)
            for key in headers[1:-1]:
                # Split the key into group and series
                group, series = key.split('_')
                # Check if the group and series exist in the detector series data
                if group in detector_series and series in detector_series[group]:
                    # Add the data value to the row
                    row_data.append(detector_series[group][series]['data'][i])
                else:
                    # Add a None or empty string depending on the output file format
                    row_data.append(None if file_extension == '.fits' else '')
            # Add the order number to the row
            row_data.append(order_number)
            # Append the row to the output data
            output_data.append(row_data)

# Sort the output data by wavelength
output_data.sort(key=lambda x: x[0])



date_str = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')

if file_extension == '.fits':
    table_data = Table(rows=output_data, names=headers)
    table_data.write(outputfilename, format='fits', overwrite=True)

    with fits.open(outputfilename, mode='update') as hdul:
        hdul[0].header['PROCSOFT'] = SWNAME
        hdul[0].header['PROCVER'] = VERSION
        hdul[0].header['DATE'] = date_str
        hdul[0].header['AUTHOR'] = AUTHOR
        hdul[0].header['EMAIL'] = EMAIL

        hdul[1].header['TUTYP1'] = 'Spectrum.Data.SpectralAxis.Value'
        hdul[1].header['TUCD1'] = 'em.wl;obs.atmos'
        hdul[1].header['TUNIT1'] = 'Angstrom'
        hdul[1].header['TDMIN1'] = min(wavelengths)
        hdul[1].header['TDMAX1'] = max(wavelengths)

elif file_extension == '.csv':
    with open(outputfilename, 'w', encoding='utf-8', newline='') as f:
        f.write(f"#PROCSOFT = {SWNAME}\n")
        f.write(f"#PROCVER = {VERSION}\n")
        f.write(f"#DATE = {date_str}\n")
        f.write(f"#AUTHOR = {AUTHOR}\n")
        f.write(f"#EMAIL = {EMAIL}\n")
        writer = csv.writer(f, lineterminator='\n') 
        writer.writerow(headers)
        writer.writerows(output_data)

input_size = os.path.getsize(inputfilename)
output_size = os.path.getsize(outputfilename)
space_gained_mb = (input_size - output_size) / (1024 * 1024)

print(f"Disk space gained: {space_gained_mb:.2f} Mb")
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time:.2f} seconds")
print("Conversion done")
print()
