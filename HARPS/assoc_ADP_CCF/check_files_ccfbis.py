#!/scratch/mbarbier/miniconda3/bin/python3
import os
import sys
import shutil
import tarfile
import glob
import requests
import warnings
import logging
import numpy as np
from astropy.io import fits


logging.basicConfig(filename='error.log', level=logging.ERROR)


# Filter the warning
warnings.filterwarnings('ignore', category=UserWarning, append=True)

# Location to store the downloaded files
storage_path = '/scratch/mbarbier/mylocal/machines/Desktop_media/SeagateHub/DATA/HARPS/harpstar/data'
cache_path = '/scratch/mbarbier/mylocal/machines/Desktop_media/SeagateHub/DATA/HARPS/harpstar/data/.cache'
exec_path='/home/mbarbier/mylocal/harpsrvcatalog/output'

file_names = ['adp', 'ccf', 'bis', 'gui']
base_filename = 'drs.{}.kwlist'

# Dictionary to store keywords and comments for different file types
data = {}

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
    data[file_type] = {'keywords': keywords_list, 'comments': comments_list}

# You can access the keywords and comments for each file type like this:
# data['adp']['keywords'] for ADP keywords
# data['adp']['comments'] for ADP comments


# Open the FITS file
with fits.open('crossmatch_ALL.fits') as hdulist:
    # Assuming the data is in the second HDU
    data = hdulist[1].data

# Initialize counters
k = 0
# Loop over all rows
for row in data:
    
    k += 1
    if k % 1000 ==0:
        print('*******',k)
    # Build the file name
    dp_id_raw = row['dp_id_raw']
    archive_id_spectra0 = row['archive_id_spectra'].strip()
    archive_id_ancillary0 = row['archive_id_ancillary'].strip()

    if len(archive_id_spectra0)>0:
        archive_id_spectra = archive_id_spectra0 + '.fits'
        dateobs = dp_id_raw[6:16]
        spectra_file = os.path.join(storage_path, dateobs, archive_id_spectra)       
        try:
            with fits.open(spectra_file) as hdulist:
                spectra_data = hdulist[0].data
        except FileNotFoundError:
            print(f'ADP file {spectra_file} not found.')
            continue
        except Exception as e:
            print(f'Error opening ADP file {spectra_file}: {e}')
            continue

        # Create patterns to search for "bis" and "ccf" files
        bis_pattern = os.path.join(storage_path, dateobs, dp_id_raw + "_bis*_A.fits")
        ccf_pattern = os.path.join(storage_path, dateobs, dp_id_raw + "_ccf*_A.fits")
        gui_pattern = os.path.join(storage_path, dateobs, dp_id_raw + "_INT_GUIDE.fits")

        # Check if the length of archive_id_ancillary is greater than 0 and
        # if there is at least one "bis" or one "ccf" file
        ccf_files = glob.glob(ccf_pattern)
        bis_files = glob.glob(bis_pattern)
        gui_files = glob.glob(gui_pattern)
        
        # Iterate through the maximum index of the two lists
        max_length = max(len(ccf_files), len(bis_files))
        for i in range(max_length):
            # If there is a ccf file at index i
            if i < len(ccf_files):
                # Open and process the i-th ccf file
                ccf_file = ccf_files[i]
                try:
                    with fits.open(ccf_file) as hdulist:
                        # Perform some operations with the ccf file
                        print(f'Processing ccf file: {ccf_file}')
                        # Example: ccf_data = hdulist[0].data
                except FileNotFoundError:
                    print(f'CCF file {ccf_file} not found.')
                    continue
                except Exception as e:
                    print(f"Error processing {ccf_file}: {e}")
                    continue
    
            # If there is a bis file at index i
            if i < len(bis_files):
                # Open and process the i-th bis file
                bis_file = bis_files[i]
                try:
                    with fits.open(bis_file) as hdulist:
                        # Perform some operations with the bis file
                        print(f'Processing bis file: {bis_file}')
                        # Example: bis_data = hdulist[0].data
                except FileNotFoundError:
                    print(f'BIS file {bis_file} not found.')
                    continue
                except Exception as e:
                    print(f"Error processing {bis_file}: {e}")
                    continue
    
            # If there is a int file at index i
            if i < len(gui_files):
                # Open and process the i-th bis file
                gui_file = gui_files[i]
                try:
                    with fits.open(gui_file) as hdulist:
                        # Perform some operations with the int file
                        print(f'Processing GUI file: {gui_file}')
                        # Example: gui_data = hdulist[0].data
                except FileNotFoundError:
                    print(f'GUI file {gui_file} not found.')
                    continue
                except Exception as e:
                    print(f"Error processing {gui_file}: {e}")
                    continue



