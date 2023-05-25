#!/usr/bin/env python3
# program that convert and UTC time in ISO format to
#   Chilean local time
#   Sideral Time
#   Local Sidereal Time assuming the location of la silla
#   MJD
import sys
from datetime import datetime
import pytz
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz

def timeconv(input_time):
    # Define time zones
    utc_tz = pytz.timezone('UTC')
    chile_tz = pytz.timezone('Chile/Continental')

    # Parse input time
    utc_time = datetime.strptime(input_time, "%Y-%m-%dT%H:%M:%S")

    # Make time timezone-aware
    utc_time = utc_tz.localize(utc_time)

    # Convert to Chile time
    chile_time = utc_time.astimezone(chile_tz)

    # Calculate MJD
    mjd_time = Time(utc_time).mjd

    # Define location for LST calculation (coordinates for La Silla Observatory in Chile)
    location = EarthLocation(lat=-29.2568, lon=-70.7392, height=2400)

    # Calculate LST
    time = Time(utc_time)
    lst_time = time.sidereal_time('mean', longitude=location.lon)
    # Calculate GMST
    gmst_time = time.sidereal_time('mean', 'greenwich')

    # Convert LST and GMST to sexagesimal format
    lst_sexagesimal = lst_time.to_string(sep=':', pad=True)
    gmst_sexagesimal = gmst_time.to_string(sep=':', pad=True)

    print(f"UTC : {utc_time.isoformat()}")
    print(f"CLT : {chile_time.isoformat()}")
    print(f"GMST: {gmst_sexagesimal}")
    print(f"LST : {lst_sexagesimal}")
    print(f"MJD : {mjd_time}")

time = sys.argv[1]
timeconv(time)
