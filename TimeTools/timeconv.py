#!/usr/bin/env python3
# program to convert times :
#   from  MJD                to  gregorian calendar
#   from  gregorian calendar to  MJD

import sys
from astropy.time import Time

def timeconv(input_time):
    try:
        # Try to convert MJD to ISO
        mjd_time = float(input_time)
        iso_time = Time(mjd_time, format='mjd').isot
        print(iso_time)
    except ValueError:
        # If conversion to float fails, assume it's an ISO date and convert to MJD
        iso_time = input_time
        mjd_time = Time(iso_time, format='isot').mjd
        print(mjd_time)

time = sys.argv[1]
timeconv(time)
