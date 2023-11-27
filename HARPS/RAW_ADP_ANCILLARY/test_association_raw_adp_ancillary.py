#!/usr/bin/python
import numpy as np
import pandas as pd
from astropy import table
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from pyvo.dal import tap

ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
tapobs = tap.TAPService(ESO_TAP_OBS)

# Function to execute ADQL query and return a pandas DataFrame
def execute_query(query):
    job = tapobs.submit_job(query)
    job.run()
    job.wait(phases=["COMPLETED", "ERROR"])
    result = job.fetch_result()
    return result.to_table().to_pandas()

# Step 1: Query the 'dbo.raw' table
raw_query = "SELECT * FROM dbo.raw WHERE instrument = 'HARPS'"
raw_data = execute_query(raw_query)

# Step 2: Filter science files from 'dbo.raw' table
science_files = raw_data[(raw_data['dp_tech'].str.startswith("ECHELLE")) & 
                         (raw_data['dp_type'].str.startswith("STAR"))]

# Step 3: Query the 'phase3v2.files' table
files_query = "SELECT * FROM phase3v2.files WHERE name LIKE '%HARPS%'"
files_data = execute_query(files_query)

# Step 4: Create new columns 'dpid1' and 'dateadp'
files_data['dpid1'] = files_data.apply(lambda row: row['name'][0:29] if row['name'].startswith("H") else row['name'].replace("TEST.", "")[0:29], axis=1)
files_data['dateadp'] = pd.to_datetime(files_data['archive_id'].str[4:30]).apply(lambda x: x.to_julian_date())

# Steps 5-7: Grouping and sorting operations
spec_group = files_data[files_data['category'] == 'SCIENCE.SPECTRUM']
anci_group = files_data[files_data['category'] == 'ANCILLARY.HARPSTAR']
spec_group_sorted = spec_group.sort_values(by='dateadp', ascending=False).drop_duplicates(subset='dpid1')
anci_group_sorted = anci_group.sort_values(by='dateadp', ascending=False).drop_duplicates(subset='dpid1')

# Step 8: Crossmatch 'dbo.raw' with 'spec'
temp1 = pd.merge(raw_data, spec_group_sorted, left_on='dp_id', right_on='dpid1', how='left')

# Step 9: Crossmatch 'temp1' with 'anci'
temp2 = pd.merge(temp1, anci_group_sorted, left_on='dp_id', right_on='dpid1', how='left', suffixes=('_spec', '_anci'))

# Step 10: Query the 'ivoa.obscore' table
obscore_query = "SELECT * FROM ivoa.obscore WHERE instrument_name = 'HARPS'"
obscore_data = execute_query(obscore_query)

# Step 11: Final crossmatch
final_table = pd.merge(temp2, obscore_data, left_on='archive_id_spec', right_on='dp_id', how='left')


