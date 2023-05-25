#!/usr/bin/env python3
# program that calculate the altazimutal coordinates
# of an object give it's equatrial coordinates
# and the date and time of observations
# does not provide accurate coordinates because
# don't take in to accout the atmospheri refraction, and
# aberration

import argparse
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
import numpy as np

def is_float(n):
    try:
        float_n = float(n)
    except ValueError:
        return False
    else:
        return True

def parse_coordinates(coord):
    if ':' in coord:  # if the coordinates are in sexagesimal format
        return SkyCoord(coord, unit=(u.hourangle, u.deg))
    else:  # if the coordinates are in degrees
        ra, dec = map(float, coord.split(","))
        return SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')

# Parse command line arguments
parser = argparse.ArgumentParser(description='Calculate altazimuthal coordinates of a star')
parser.add_argument('-L', '--label', default='OBJECT', help='Star name')
parser.add_argument('-c', '--coordinates', nargs=2, type=float, required=True, help='Star equatorial coordinates')
parser.add_argument('-s', '--site', default='-29.2568,-70.7392,2400', help='Site coordinates (latitude,longitude,altitude in meters)')
parser.add_argument('-t', '--time', required=True, help='Time (JD, MJD, or ISO format)')

args = parser.parse_args()

# Parse coordinates
ra, dec = args.coordinates

# Parse site coordinates
observer_lat, observer_lon, observer_alt = map(float, args.site.split(","))

#Parse time
if is_float(args.time):
    time_float = float(args.time)
    if time_float > 1000000:
        # Treat as Julian Date
        time = Time(time_float, format='jd')
    else:
        # Treat as Modified Julian Date
        time = Time(time_float, format='mjd')
else:
    # Treat as ISO formatted string
    if "T" in args.time:
        time = Time(args.time, format='isot')
    else:
        time = Time(args.time, format='iso')

# Define star and observer location
star = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')

# Define observer location
observer_location = EarthLocation(lat=observer_lat*u.deg, lon=observer_lon*u.deg, height=observer_alt*u.m)

# Transform to AltAz
altaz = star.transform_to(AltAz(obstime=time, location=observer_location))

# Print altazimuthal coordinates
print(f"Altazimuthal coordinates of {args.label}:")
print(f"Altitude : {altaz.alt.deg} degrees")
print(f"Azimuth  : {altaz.az.deg} degrees")


# Check if the star is rising or setting
if altaz.az.deg < 180:
    print(f"{args.label} is rising.")
else:
    print(f"{args.label} is setting.")
