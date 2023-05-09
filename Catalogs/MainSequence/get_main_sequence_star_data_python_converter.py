#!/usr/bin/env python3
import pandas as pd
from scipy.interpolate import interp1d
from scipy.interpolate import Rbf
import sys

# Read CSV file into a pandas DataFrame
input_file = 'msdata.csv'
output_file = 'main_sequence_data.csv'
df = pd.read_csv(input_file)
df.columns = df.columns.str.strip()

#sys.exit()

# Identify the rows where 'logL' is missing
# Get the rows where 'logL' is not missing
# Create a interpolation/extrapolation function using the rows without missing 'logL'
# Fill missing values in the 'logL' column using the interpolation/extrapolation function
missing_values_mask_logL = df['logL'].isna()
not_missing_values_mask_logL = ~missing_values_mask_logL
f_logL = interp1d(df.loc[not_missing_values_mask_logL, 'Teff'], df.loc[not_missing_values_mask_logL, 'logL'], kind='linear', fill_value='extrapolate')
df.loc[missing_values_mask_logL, 'logL'] = f_logL(df.loc[missing_values_mask_logL, 'Teff'])

# Identify the rows where 'Msun' is missing
# Get the rows where 'Msun' is not missing
# Create a interpolation/extrapolation function using the rows without missing 'Msun'
# Fill missing values in the 'Msun' column using the interpolation/extrapolation function
missing_values_mask_Msun = df['Msun'].isna()
not_missing_values_mask_Msun = ~missing_values_mask_Msun
f_Msun = interp1d(df.loc[not_missing_values_mask_Msun, 'Teff'], df.loc[not_missing_values_mask_Msun, 'Msun'], kind='linear', fill_value='extrapolate')
df.loc[missing_values_mask_Msun, 'Msun'] = f_Msun(df.loc[missing_values_mask_Msun, 'Teff'])


# Iterate through each column (excluding the first and second columns)
for column in df.columns[2:]:
	# Drop rows with missing values in the current column and Teff column
	temp_df = df[['Teff', column]].dropna()

	# Perform linear interpolation/extrapolation using Teff as the X variable
	f = interp1d(temp_df['Teff'], temp_df[column], kind='linear', fill_value='extrapolate')

	# Fill missing values in the column using the extrapolation function
	missing_values_mask = df[column].isna()
	df.loc[missing_values_mask, column] = f(df.loc[missing_values_mask, 'Teff'])

# Write the updated DataFrame to a new CSV file
df.to_csv(output_file, index=False)


sys.exit()



# seems nice but in reality don't work
if (1 == 2):
    # Define the independent variables
    independent_vars = df[['Teff', 'logL', 'Msun']]
    # Iterate over the other columns from column 4 to the one before the last, excluding column 5
    for column in df.columns[4:-1]:
        if column != 'logL':
            # Create mask for non-missing values in the current column
            non_missing_values_mask = df[column].notna()
            
            # Create the RBF interpolator using the independent variables and the current column
            rbf = Rbf(independent_vars.loc[non_missing_values_mask, 'Teff'], 
                      independent_vars.loc[non_missing_values_mask, 'logL'], 
                      independent_vars.loc[non_missing_values_mask, 'Msun'], 
                      df.loc[non_missing_values_mask, column], 
                      function='linear')
    
            # Create mask for missing values in the current column
            missing_values_mask = df[column].isna()
    
            # Use the RBF interpolator to fill missing values in the current column
            df.loc[missing_values_mask, column] = rbf(independent_vars.loc[missing_values_mask, 'Teff'], 
                                                       independent_vars.loc[missing_values_mask, 'logL'], 
                                                       independent_vars.loc[missing_values_mask, 'Msun'])






