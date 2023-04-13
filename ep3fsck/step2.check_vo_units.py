#!/bin/python3
import pandas as pd
import pint
from pint.errors import UndefinedUnitError

ureg = pint.UnitRegistry()
ureg.define('angstrom = 1e-10 * meter = \u00c5')

# Define a function to check the units
def check_units(unit_str):
    try:
        ureg(unit_str)
    except UndefinedUnitError:
        return f"Invalid unit: {unit_str}"
    except:
        return f"Error checking unit: {unit_str}"

# Read the CSV table into a pandas DataFrame
df = pd.read_csv('keywords_vo_table.csv')

# Apply the check_units function to the TUNIT column and add a new column for the error message
df['unit_error'] = df['TUNIT'].apply(check_units)

# Filter the DataFrame to show only the rows with unit errors
unit_errors = df.loc[df['unit_error'] != '', ['COLUMN', 'TUNIT', 'unit_error']]

# Print the unit errors
print(unit_errors)
