#!/bin/bash
url="https://dataportal.eso.org/dataportal_new/file/"
wdir=$PWD
mylist="path_and_filename_of_DPID_to_retrieve"  
# This file need to contain only 1 column
# the content of the column need to be a valid DPID of a ESO archive file

function getfits() {
# Count the total number of lines in the mylist file
N=$(wc -l < "$mylist")

# Initialize the counter and the last progress variable
i=0
last_progress=-10

# Read the mylist file line by line
while read -r dp_id; do
# Increment the counter
((i++))

# Calculate the percentage of progress
percentage=$((i * 100 / N))

# Check if the progress is a multiple of 10%
if ((percentage / 10 > last_progress / 10)); then
echo "Progress: $percentage%"
last_progress=$percentage
fi

# Create the directory if it doesn't exist
outdir="telescope/instrument/data_type"  # adjust the path
mkdir -p "$outdir"

# Download the file
if ! wget -q -O "$outdir/$dp_id.fits" "$url$dp_id"; then
echo "Failed to download $dp_id"
continue
fi

# Uncompress the file if it's compressed
file_type=$(file -b "$outdir/$dp_id.fits")
if [[ $file_type == "compress"* ]]; then
cd "$outdir"
mv "$dp_id.fits" "$dp_id.fits.Z"
uncompress "$dp_id.fits.Z"
cd "$wdir"
fi

# Create a header file for the FITS file
dfits -x 0 "$outdir/$dp_id.fits" > "$outdir/$dp_id.hdr"

done < "mylist"

}

getfits

