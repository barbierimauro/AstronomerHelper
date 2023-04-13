#!/bin/python3
from astropy.io.votable.ucd import parse_ucd, check_ucd
import pandas as pd
# Load the CSV table into a pandas DataFrame
df = pd.read_csv('keywords_vo_table.csv')

def check_and_suggest(idx, ucd):
	try:
		parse_ucd(ucd,check_controlled_vocabulary=True)
		print(f"COLUMN {idx+1}, UCD: '{ucd}', CORRECT")
	except ValueError as err:
	# Print the error message and the corresponding COLUMN value
		print(f"COLUMN {idx+1}, UCD: '{ucd}', INVALID, error message: {err}")

# Loop through the TUCD column and check the syntax of each UCD
for idx, ucd in df['TUCD'].items():
	check_and_suggest(idx, ucd)
