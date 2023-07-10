#!/bin/bash
# 2023 maurobarbieri.science@gmail.com
#
# This script retrieves the headers of an Advanced Data Product (ADP) file 
# given its Data Product Identifier (DPID). It then scans the headers for all 
# occurrences of the PROVENANCE keyword and retrieves the corresponding headers 
# for each provenance. Lastly, it searches the provenance headers for user-defined 
# keywords or patterns.
#
# Usage: 
#   hdr_info.sh DPID_ADP
#
# Arguments:
#   DPID_ADP: The Data Product Identifier of an ADP file.
#
# Dependencies:
#   This script requires the following commands to be available on your system:
#     - wget
#     - sed
#     - grep
#     - cut
#     - tr
#
# User-defined Settings:
#   Modify the "keywords_to_search" variable to add more keywords or other 
#   patterns to search in the headers. Keywords should be separated by the 
#   pipe character ("|").
#
# Output:
#   For each provenance found, this script outputs a header file in the format:
#       {provenance}.{dpid_provenance}.hdr
#   It also prints out lines from the header that match any of the user-defined 
#   keywords or patterns.
#




# ADD YOUR KEYWORD TO SEARCH SEPARATED BY "|"
keywords_to_search='OBJECT  =|DATE-OBS='


# Exit script on any error
set -e

# url headers
url_hdr="http://archive.eso.org/hdr?DpId="

# Check if the script received an argument
if [ -z "$1" ]; then
    echo "Usage: $0 DPID_ADP"
    exit 1
fi

dpid_adp=$1
# convert ":" to "_" in filenames
dpid_adp_nocolons=${dpid_adp//:/_}
adp=$url_hdr$dpid_adp

# retrieve the ADP header
wget -q -O "${dpid_adp_nocolons}.tmp" $adp > /dev/null
# remove the HTML part
sed -n '/<pre>/,/<\/pre>/p' "${dpid_adp_nocolons}.tmp" | sed '1d;$d' > "${dpid_adp_nocolons}.hdr"
rm "${dpid_adp_nocolons}.tmp"
# print basic information
echo
echo "ADP input file"
echo $dpid_adp   "    :   " $dpid_adp
echo "header in : " "${dpid_adp}.hdr"
echo
echo "*___________________________________________________________*"
# start the loop on the PROVENANCE keywords
# grep for PROV and read the line
grep "PROV" "${dpid_adp_nocolons}.hdr" | while read -r line_adp
do
    # get the keyword name
    provenance=$(echo $line_adp | cut -d= -f1 | tr -d ' ')
    # get the keyword content
    dpid_provenance=$(echo $line_adp | cut -d= -f2 | cut -d/ -f1 | tr -d ' '| tr -d "'" | sed -e "s/.fits//g")
    raw=$url_hdr$dpid_provenance
    echo $provenance "    :   " $dpid_provenance
    # get the PROVENANCE file
    wget -q -O "${provenance}.${dpid_provenance}.tmp" $raw > /dev/null
    # remove the HTML part
    sed -n '/<pre>/,/<\/pre>/p' "${provenance}.${dpid_provenance}.tmp" | sed '1d;$d' > "${provenance}.${dpid_provenance}.hdr"
    rm "${provenance}.${dpid_provenance}.tmp"

    echo "header in : " "${provenance}.${dpid_provenance}.hdr"
    echo
    # search for the relevant keywords and/or patterns
    grep -hiE "$keywords_to_search" "${provenance}.${dpid_provenance}.hdr" | while read -r line_raw
    do
       # print the content
       echo $line_raw
    done
    echo "*___________________________________________________________*"
    echo
done


