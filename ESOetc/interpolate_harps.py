#!/usr/bin/python3
import pandas as pd
import numpy as np

# Load the data
df = pd.read_csv('harps.csv')

# Get the range of wavelengths
min_wl = df['wl'].min()
max_wl = df['wl'].max()

# Create new wavelength scale
new_wl = np.arange(min_wl, max_wl, 0.01)
new_wl = np.around(new_wl, 2)

# Initialize new DataFrame with the new wavelength scale
df_new = pd.DataFrame(new_wl, columns=['wl'])

# Interpolate other columns
for col in df.columns[1:]:
    df_new[col] = np.interp(new_wl, df['wl'], df[col])

# Save the new DataFrame to a new CSV file
df_new.to_csv('harps_out.csv', index=False)
