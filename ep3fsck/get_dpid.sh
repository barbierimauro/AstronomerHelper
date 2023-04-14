#!/bin/bash
#
# Script for downloading multiple files from the ESO archive
# It is necessary to known the name of the files and have stored them
# in a list one per line


#BEGIN PERSONAL SETUP

# modify this variable with your file name
mylist="path_and_filename_of_DPID_to_retrieve"  
# This file need to contain only 1 column
# the content of the column need to be a valid DPID of a ESO archive file

# adjust the path
outdir="telescope/instrument/data_type"  

# you can modify this number in integers from 1 to 99
# after downloading this precentage of files over the total number 
# it will be displayed a notification on the screen
# no need to modify in general case, but useful when downloading thousands of files
percento=10

if ! [[ $percento =~ ^[0-9]+$ ]]; then
    echo "Error: percentual is not numeric"
elif (( percento < 1 || percento > 99 )); then
    echo "Error: percentual is not between 1 and 99"
fi


#END PERSONAL SETUP

url="https://dataportal.eso.org/dataportal_new/file/"
wdir=$PWD

function getfits() {
# Count the total number of lines in the mylist file
N=$(wc -l < "$mylist")

# Initialize the counter and the last progress variable
i=0
last_progress=-$percento

# Read the mylist file line by line
while read -r dp_id; do
# Increment the counter
((i++))

# Calculate the percentage of progress
percentage=$((i * 100 / N))

# Check if the progress is a multiple of $percento 
if ((percentage / $percento > last_progress / $percento)); then
echo "Progress: $percentage%"
last_progress=$percentage
fi

# Create the directory if it doesn't exist
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
# not needed in normal operations, but can be activated uncommenting the following line
# dfits -x 0 "$outdir/$dp_id.fits" > "$outdir/$dp_id.hdr"

done < "mylist"

}

getfits

