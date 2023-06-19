#!/usr/bin/python3
from astropy.io import fits
import csv

# Define the input and output file paths
input_file = 'tab_ccfbis.csv'
output_file = 'tab_ccfbis.output.csv'
adp_kwlist = 'drs.adp.kwlist'
ccf_kwlist = 'drs.ccf.kwlist'
bis_kwlist = 'drs.bis.kwlist'

# Read the list of FITS keywords from the txt files
with open(adp_kwlist, 'r') as f:
    adp_keywords_list = [line.strip() for line in f]

with open(ccf_kwlist, 'r') as f:
    ccf_keywords_list = [line.strip() for line in f]

with open(bis_kwlist, 'r') as f:
    bis_keywords_list = [line.strip() for line in f]

# Create the output CSV file and write the header row
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['rawfile', 'dateobs', 'mjdobs', 'yearobs', 'mask', 'flag_cb'] + ccf_keywords_list + bis_keywords_list)

    # Loop over each line in the input CSV file
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Extract the values of DAT1
            dat1 = [row['rawfile'], row['dateobs'], row['mjdobs'], row['yearobs'], row['mask'], row['flag_cb']]

            # Extract the values of DAT2 from the ccffile
            with fits.open(row['ccffile']) as hdul:
                header = hdul[0].header
                dat2 = [header.get(keyword, '') for keyword in ccf_keywords_list]

            # Extract the values of DAT3 from the bisfile
            with fits.open(row['bisfile']) as hdul:
                header = hdul[0].header
                dat3 = [header.get(keyword, '') for keyword in bis_keywords_list]

            # Write the row to the output CSV file
            writer.writerow(dat1 + dat2 + dat3)



