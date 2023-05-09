#!/usr/bin/python3 
'''
Title: 
Stellar Catalog Generator for Testing Cross Correlation Algorithms

Purpose: 
The purpose of this program is to create test catalogs for 
performing catalog cross-correlation tests. It generates a synthetic 
catalog of reference stars and their companion stars with desired 
properties. These test catalogs can be used to compare and analyze 
various catalog matching and cross-correlation techniques.

Program Overview:
The program generates two catalogs:

-A reference catalog (catalog_ref), containing the main set of reference stars with their properties.
-A comparison catalog (catalog_cmp), containing three types of stars:
   - The reference stars with added variability in position and parameters values.
   - Companion stars generated within a specified radius around the reference stars.
   - Background stars, which are unrelated to the reference stars but generated within the same region.

Inputs:

    N: Number of reference stars.
    params_ref_cat: A list of tuples, where each tuple contains the mean, standard deviation, and error for each variable in the reference catalog.
    params_2nd_cat: A list of tuples, where each tuple contains the minimum and maximum values for each variable in the comparison catalog.
    radius: The radius within which companion stars are generated (in arcsec).
    min_M: Minimum number of background stars to be generated.
    max_M: Maximum number of background stars to be generated.
    lon_sig: Longitude error values in mas.
    lat_sig: Latitude error values in mas.
    max_var_percen: Maximum percentage of variability allowed for each variable.
    proba: Probability distribution used to generate the number of companions.
    min_M: Minimum number of background stars to be generated.
    max_M: Maximum number of background stars to be generated.
    params_cmp: A list of tuples defining the parameters for the comparison catalog.
    ncomp: Number of companion stars.
    seed: Seed for the random number generator (optional).

Outputs:

    catalog_ref: A pandas DataFrame containing the reference catalog with the following columns:
        id_ref: Reference star ID.
        ra: Right ascension (in degrees).
        dec: Declination (in degrees).
        e_ra: Right ascension error (in mas).
        e_dec: Declination error (in mas).
        plx, pmra, pmde, rv, vmag, vi: Variables for each star with their respective errors 
        (e_plx, e_pmra, e_pmde, e_rv, e_vmag, e_vi).
        
    catalog_cmp: A pandas DataFrame containing the comparison catalog with the following columns:
        id_cmp: Comparison star ID.
        id_cmp_ref: ID of the reference star associated with the comparison star.
        id_cmp_label: Label for the type of comparison star 
        (A for reference stars with variability, B-Z for companion stars, or 'x' for background stars).
        ra: Right ascension (in degrees).
        dec: Declination (in degrees).
        e_ra: Right ascension error (in mas).
        e_dec: Declination error (in mas).
        plx, pmra, pmde, rv, vmag, vi: Variables for each star with their respective errors 
        (e_plx, e_pmra, e_pmde, e_rv, e_vmag, e_vi).
        
'''

import numpy as np
import pandas as pd
import time
import sys
from astropy.table import Table
from astropy.table import QTable
from astropy.table import Column
from astropy.io import fits


# Function to generate random points on a sphere, with error values for longitude and latitude
def random_point_on_sphere(N, lon_sig, lat_sig):
    lon = 2 * np.pi * np.random.rand(N)  # Random longitude values
    sin_lat = 2 * np.random.rand(N) - 1  # Random sine latitude values
    lat = np.arcsin(sin_lat)  # Calculate latitude values from sine latitude
    lon_err = np.abs(np.random.normal(0, lon_sig, N))  # Longitude error values
    lat_err = np.abs(np.random.normal(0, lat_sig, N))  # Latitude error values
    coord = np.array([lon, lat, lon_err, lat_err]).T  # Combine coordinates and errors into an array
    return coord

# Function to generate synthetic catalogs of reference and companion stars
def generate_ref_data(N, lon_sig, lat_sig, params, radius, max_var_percen, proba, min_M, max_M, params_cmp, ncomp, seed=None):
    # Set the random seed if provided, otherwise use the current time
    if seed is None:
        seed = int(time.time())
    np.random.seed(seed)
    
    num_vars = len(params)  # Number of variables for each star
    values = np.zeros((N, num_vars))  # Array to store variable values for reference stars
    errors = np.zeros((N, num_vars))  # Array to store variable errors for reference stars
    
    # Define the columns and data types for the reference and companion catalogs
    catalog_ref_columns = ['id_ref', 'ra', 'de', 'e_ra', 'e_de', 'plx', 'pmra', 'pmde', 'rv', 'vmag', 'vicol', 'e_plx', 'e_pmra', 'e_pmde', 'e_rv', 'e_vmag', 'e_vicol']
    catalog_cmp_columns = ['id_cmp', 'id_cmp_ref', 'id_cmp_label', 'ra', 'de', 'e_ra', 'e_de', 'plx', 'pmra', 'pmde', 'rv', 'vmag', 'vicol', 'e_plx', 'e_pmra', 'e_pmde', 'e_rv', 'e_vmag', 'e_vicol']
    
    catalog_ref = pd.DataFrame(columns=catalog_ref_columns)  # Create an empty DataFrame for the reference catalog
    catalog_cmp = pd.DataFrame(columns=catalog_cmp_columns)  # Create an empty DataFrame for the companion catalog


     # Dictionary to store units of each column in the reference catalog
    units_ref = {'id_ref': 'int', 
             'ra': 'deg', 'de': 'deg', 'e_ra': 'mas', 'e_de': 'mas',
             'plx': 'mas', 'pmra': 'mas/yr', 'pmde': 'mas/yr', 'rv': 'km/s', 'vmag': 'mag', 'vicol': 'mag', 
             'e_plx': 'mas', 'e_pmra': 'mas/yr', 'e_pmde': 'mas/yr', 'e_rv': 'km/s', 'e_vmag': 'mag', 'e_vicol': 'mag'}

    # Dictionary to store units of each column in the comparison catalog
    units_cmp = {'id_ref': 'int', 'id_cmp_ref': 'int', 'id_cmp_label': 'str', 
             'ra': 'deg', 'de': 'deg', 'e_ra': 'mas', 'e_de': 'mas',
             'plx': 'mas', 'pmra': 'mas/yr', 'pmde': 'mas/yr', 'rv': 'km/s', 'vmag': 'mag', 'vicol': 'mag', 
             'e_plx': 'mas', 'e_pmra': 'mas/yr', 'e_pmde': 'mas/yr', 'e_rv': 'km/s', 'e_vmag': 'mag', 'e_vicol': 'mag'}
    
    # Dictionary to store data types of each column in the reference catalog
    data_types_ref = {'id_ref': int, 
             'ra': float, 'de': float, 'e_ra': float, 'e_de': float,
             'plx': float, 'pmra': float, 'pmde': float, 'rv': float, 'vmag': float, 'vicol': float, 
             'e_plx': float, 'e_pmra': float, 'e_pmde': float, 'e_rv': float, 'e_vmag': float, 'e_vicol': float}
    
    # Dictionary to store data types of each column in the comparison catalog
    data_types_cmp = {'id_cmp': int, 'id_cmp_ref': int, 'id_cmp_label': str, 
             'ra': float, 'de': float, 'e_ra': float, 'e_de': float,
             'plx': float, 'pmra': float, 'pmde': float, 'rv': float, 'vmag': float, 'vicol': float, 
             'e_plx': float, 'e_pmra': float, 'e_pmde': float, 'e_rv': float, 'e_vmag': float, 'e_vicol': float}

    # Assign the units dictionaries to the catalog DataFrames
    catalog_ref.attrs['units'] = units_ref
    catalog_cmp.attrs['units'] = units_cmp
    
    # Set the data types for each catalog DataFrame
    catalog_ref = catalog_ref.astype(data_types_ref)
    catalog_cmp = catalog_cmp.astype(data_types_cmp)

    #print(catalog_ref.dtypes)
    #print(catalog_cmp.dtypes)


    # Print a message indicating the generation of reference stars
    print('reference stars')
    
    # Loop through the reference stars
    for n in range(N):
        # Generate random coordinates on a sphere
        coord = random_point_on_sphere(1, lon_sig, lat_sig)
        
        # Assign an ID to the reference star
        id_ref = int(n + 1)
        
        # Generate coordinates with added variability
        coord2 = coord.copy() * (1 + np.random.uniform(-max_var_percen, max_var_percen, coord.shape))
        
        # Loop through the variables
        for i in range(num_vars):
            mu, sigma, err = params[i]
            values[n, i] = np.random.normal(mu, sigma)
            errors[n, i] = np.random.normal(0, err)
            
            # Generate values and errors with added variability
            values2 = values.copy() * (1 + np.random.uniform(-max_var_percen, max_var_percen, values.shape))
            errors2 = errors.copy() * (1 + np.random.uniform(-max_var_percen, max_var_percen, errors.shape))
        
        # Add the generated data to the reference and comparison catalogs
        catalog_ref.loc[n] = np.hstack([int(id_ref), coord.ravel(), values[n, :], errors[n, :]])
        catalog_cmp.loc[n] = np.hstack([int(id_ref), int(id_ref), 'A', coord2.ravel(), values2[n, :], errors2[n, :]])
    
    # Print a message indicating the generation of companion stars
    print('companion stars')
    
    # Initialize the index for companion stars
    n1 = n + 1
    
    # Loop through the reference stars
    for n in range(N):
        # Generate a random number of companions for each reference star
        rand_num = np.random.uniform(0, 1)
        L = ncomp[np.searchsorted(proba, rand_num)]
        
        # Get the coordinates of the reference star
        lon_center = catalog_ref.iloc[n, 1]
        lat_center = catalog_ref.iloc[n, 2]

        # Generate companion stars if there are any
        if (L > 0):
            # Initialize arrays for the values and errors of the companion stars
            values_2nd = np.zeros((L, num_vars))
            errors_2nd = np.zeros((L, num_vars))
    
            # Loop through the companion stars
            for i in range(L):
                n1=n1+1
                theta = 2 * np.pi * np.random.rand()
                rho = np.abs(np.random.normal(0, radius))
                lat = np.arcsin(np.sin(lat_center) * np.cos(rho) + np.cos(lat_center) * np.sin(rho) * np.cos(theta))
                lon = lon_center + np.arctan2(np.sin(theta) * np.sin(rho) * np.cos(lat_center), np.cos(rho) - np.sin(lat_center) * np.sin(lat))
                lon_err = np.abs(np.random.normal(0, lon_sig))
                lat_err = np.abs(np.random.normal(0, lat_sig))

                # Generate the label for the companion star
                label = chr(66 + i) if i <= 25 else 'Z'

                # Create the coordinates for the companion star
                coord_2nd = np.array([lon, lat, lon_err, lat_err])

                # # Loop through the variables
                # for j in range(num_vars):
                    # # Generate values and errors with added variability
                    # mu, sigma, err = params[j]
                    # values_2nd[i, j] = np.random.normal(mu, sigma)
                    # errors_2nd[i, j] = np.random.normal(0, err)

                # Generate the companion star's variable values and errors
                # Loop through the variables
                random_numbers=np.random.uniform(size=6)
                for j in range(num_vars):
                    # Initialize the values with NaN
                    values_2nd[i, j] = np.nan
                    errors_2nd[i, j] = np.nan
                    
                    if j == 0:
                        random_number = random_numbers[j]
#                        print(j,random_number)
                        if random_number > 0.1:
                            # Generate first value : plx
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
                    elif j == 1:
                        random_number = random_numbers[0]
#                        print(j,random_number)
                        if random_number > 0.1:
                            # Generate second value : pmra
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
                    elif j == 2:
                        random_number = random_numbers[0]
#                        print(j,random_number)
                        if random_number > 0.1:
                            # Generate third value : pmde
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
                    elif j == 3: 
                        random_number = random_numbers[j]
#                        print(j,random_number)
                        if random_number > 0.7:
                            # Generate the fourth value : vrad
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
                    elif j == 4 :
                        random_number = random_numbers[j]
                        #print(j,random_number)
                        if random_number > 0.1:
                            # Generate the fifth value : vmag
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
                            if values_2nd[i, j]>15:
                                 values_2nd[i, 3] = np.nan
                                 errors_2nd[i, 3] = np.nan
                    elif j == 5 and not np.isnan(values_2nd[i, 4]):
                        random_number = random_numbers[j]
#                        print(j,random_number)
                        if random_number > 0.2:
                            # Generate the sixth value : v-i
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
#                print(n,i,values_2nd[i, 0:6])            

            # Add the generated data to the comparison catalog
            catalog_cmp.loc[n1] = np.hstack([int(n + 1), int(n1 + 1), label, coord_2nd.ravel(), values_2nd[i, :], errors_2nd[i, :]])
            #print(label,type(label),catalog_cmp['id_cmp_label'].dtype)
            
# Print statement to indicate the generation of background stars
    print('background stars')
    # Get the number of records in the catalog_cmp DataFrame
    n2= len(catalog_cmp)
    # Iterate through all reference stars
    for n in range(N):
    # Generate a random number of background stars M, within the range of min_M and max_M
        M = np.random.randint(min_M, max_M + 1)

    # Get the reference star's longitude and latitude
        lon_center = catalog_ref.iloc[n, 1]
        lat_center = catalog_ref.iloc[n, 2]

    # Check if there are background stars to be generated
        if (M > 0):
        # Initialize arrays to store generated values and errors for background stars
            values_2nd = np.zeros((M, num_vars))
            errors_2nd = np.zeros((M, num_vars))

        # Generate background stars
            for i in range(M):
                n2=n2+1
                # Generate random polar coordinates for the background star
                theta = 2 * np.pi * np.random.rand()
                rho = radius*np.random.rand()
                # Convert polar coordinates to geographic coordinates (longitude and latitude)
                lat = np.arcsin(np.sin(lat_center) * np.cos(rho) + np.cos(lat_center) * np.sin(rho) * np.cos(theta))
                lon = lon_center + np.arctan2(np.sin(theta) * np.sin(rho) * np.cos(lat_center), np.cos(rho) - np.sin(lat_center) * np.sin(lat))

            # Generate longitude and latitude errors using the given standard deviations
                lon_err = np.abs(np.random.normal(0, lon_sig))
                lat_err = np.abs(np.random.normal(0, lat_sig))

            # Generate the label for the background star
                label = 'x'

            # Create an array with the background star's coordinates and errors
                coord_2nd = np.array([lon, lat, lon_err, lat_err])

            # # Generate the background star's variable values and errors
                # for j in range(num_vars):
                    # vmin, vmax = params_cmp[j]
                    # err=np.abs(vmax-vmin)*0.02
                    # values_2nd[i, j] = np.random.uniform(vmin, vmax)
                    # errors_2nd[i, j] = np.random.uniform(0,err)

            # Generate the background star's variable values and errors
                # Loop through the variables
                random_numbers=np.random.uniform(size=6)
                for j in range(num_vars):
                    # Initialize the values with NaN
                    values_2nd[i, j] = np.nan
                    errors_2nd[i, j] = np.nan
                    
                    if j == 0:
                        random_number = random_numbers[j]
#                        print(j,random_number)
                        if random_number > 0.1:
                            # Generate first value : plx
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
                    elif j == 1:
                        random_number = random_numbers[0]
#                        print(j,random_number)
                        if random_number > 0.1:
                            # Generate second value : pmra
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
                    elif j == 2:
                        random_number = random_numbers[0]
#                        print(j,random_number)
                        if random_number > 0.1:
                            # Generate third value : pmde
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
                    elif j == 3: 
                        random_number = random_numbers[j]
#                        print(j,random_number)
                        if random_number > 0.7:
                            # Generate the fourth value : vrad
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
                    elif j == 4 :
                        random_number = random_numbers[j]
                        #print(j,random_number)
                        if random_number > 0.1:
                            # Generate the fifth value : vmag
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
                            if values_2nd[i, j]>15:
                                 values_2nd[i, 3] = np.nan
                                 errors_2nd[i, 3] = np.nan
                    elif j == 5 and not np.isnan(values_2nd[i, 4]):
                        random_number = random_numbers[j]
#                        print(j,random_number)
                        if random_number > 0.2:
                            # Generate the sixth value : v-i
                            mu, sigma, err = params[j]
                            values_2nd[i, j] = np.random.normal(mu, sigma)
                            errors_2nd[i, j] = np.random.normal(0, err)
#                print(n,i,values_2nd[i, 0:6])            

                # Add the background star to the catalog_cmp DataFrame
                catalog_cmp.loc[n2] = np.hstack([int(n2), int(n+1) , label, coord_2nd.ravel(), values_2nd[i, :], errors_2nd[i, :]])
                #print(label,type(label),catalog_cmp['id_cmp_label'].dtype)
                
    radian_to_deg = 180.0/np.pi
#    radian_to_arcsec = 180.0*3600.0/np.pi
    
    # Convert the reference and companion catalogs' longitude and latitude from radians to degrees

    catalog_ref['ra']  *= radian_to_deg
    catalog_ref['de'] *= radian_to_deg
    catalog_cmp['ra'] = catalog_cmp['ra'].astype(float) * radian_to_deg
    catalog_cmp['de'] = catalog_cmp['de'].astype(float) * radian_to_deg

    catalog_ref = catalog_ref.astype(data_types_ref)
    catalog_cmp = catalog_cmp.astype(data_types_cmp)

     # Return the generated reference and companion catalogs

    return catalog_ref,catalog_cmp


# parameters    

radian_to_arcsec = 180.0*3600.0/np.pi
radian_to_deg = 180.0/np.pi


# Define the number of reference stars (N)
N = 100

# Define the parameters for the reference catalog (params_ref_cat)
# Each tuple contains the mean, standard deviation, and error for each variable
#                  plx          pmra         pmde         rv          vmag         vi
params_ref_cat = [(50, 50, 5), (0, 200, 5), (0, 200, 5), (0, 50, 1), (10, 5, 0.1), (1.5, 1.0, 0.1)]

# Define the parameters for the comparison catalog (params_2nd_cat)
# Each tuple contains the minimum and maximum values for each variable
params_2nd_cat = [(-10, 700), (-4000, 4000), (-4000, 4000), (-100, 100), (0, 22), (-1, 6)]

# Define the radius within which companion stars are generated (radius)
# Convert the radius from arcseconds to radians
radius = 10 / radian_to_arcsec

# Define the minimum and maximum number of background stars to be generated (min_M, max_M)
min_M = 1
max_M = 20

# Define the longitude and latitude error values (lon_sig, lat_sig)
# Convert the error values from mas (milliarcseconds) to radians
lon_sig = 200.0 #/ (radian_to_arcsec * 1000)
lat_sig = 200.0 #/ (radian_to_arcsec * 1000)

# Define the maximum percentage of variability allowed for each variable (max_var_percen)
max_var_percen = 1e-5

# Define the probability distribution for the number of companions (proba)
# and the corresponding number of additional companions (ncomp)
proba=[0.56,0.89,0.97,1.00]
ncomp=[0   ,1   ,2   ,3   ] 
# number of additonal companion to create
# 1 = binary     star (1+1)
# 2 = ternary    star (1+2)
# 1 = quaternary star (1+3)

# Define the random seed for reproducibility (seed)
#seed = 69
seed = None

# Call the function to generate reference and companion catalogs
catalog_ref,catalog_cmp = generate_ref_data(N, lon_sig, lat_sig, params_ref_cat, radius,max_var_percen, proba, min_M, max_M, params_2nd_cat, ncomp, seed)

# Save the catalogs to CSV files
catalog_ref.astype({'id_ref': int}).to_csv('catalog_ref.csv', index=False)
catalog_cmp.to_csv('catalog_cmp.csv', index=False)

# Convert DataFrames to Astropy Tables
#ref_table = Table.from_pandas(catalog_ref.astype({'id_ref': int}))
#cmp_table = Table.from_pandas(catalog_cmp)
#ref_table = QTable.from_pandas(catalog_ref.astype({'id_ref': int}))
#cmp_table = QTable.from_pandas(catalog_cmp)

# Add the columns one by one with the desired data type
ref_table = QTable([Column(catalog_ref[col].astype(dtype), name=col) for col, dtype in zip(catalog_ref.columns, [int, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float])])
cmp_table = QTable([Column(catalog_cmp[col].astype(dtype), name=col) for col, dtype in zip(catalog_cmp.columns, [int, int, str, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float])])


#print(catalog_ref.dtypes)
#print(catalog_cmp.dtypes)
#sys.exit()

# Define UCD and units keywords for the columns
ucd_keywords_ref = ['ID_MAIN', 'POS_EQ_RA_MAIN', 'POS_EQ_DEC_MAIN', 'PHOT_MAG_R', 'PHOT_COLOR_G-R', 'PHYS_SIZE_RADIUS_GENERAL']
units_keywords_ref = [None, 'deg', 'deg', 'mag', 'mag', 'arcsec']
ucd_keywords_cmp = ['ID_MAIN', 'POS_EQ_RA_MAIN', 'POS_EQ_DEC_MAIN', 'PHOT_MAG_R', 'PHOT_COLOR_G-R', 'PHYS_SIZE_RADIUS_GENERAL']
units_keywords_cmp = [None, 'deg', 'deg', 'mag', 'mag', 'arcsec']

# Add UCD and units keywords to the headers
for i, (ucd, unit) in enumerate(zip(ucd_keywords_ref, units_keywords_ref)):
    ref_table.columns[i].meta = {'UCD': ucd, 'unit': unit}

for i, (ucd, unit) in enumerate(zip(ucd_keywords_cmp, units_keywords_cmp)):
    cmp_table.columns[i].meta = {'UCD': ucd, 'unit': unit}

sys.exit()
# Save the catalogs as FITS binary tables
ref_table.write('catalog_ref.fits', format='fits', overwrite=True)
cmp_table.write('catalog_cmp.fits', format='fits', overwrite=True)

