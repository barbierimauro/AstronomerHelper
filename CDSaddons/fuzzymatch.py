#!/bin/python3
import re
import unicodedata
from fuzzywuzzy import process

def remove_special_characters(text):
    chars_to_remove = r'[]!~@#${}()*^%_=;\'`?/\\'
    return re.sub(f'[{re.escape(chars_to_remove)}]', '', text)

def find_closest_catalog(catalogs, text):
    return process.extractOne(text, catalogs)[0]

def remove_latex_symbols(text):
    return re.sub(r'\\(?:text)?[a-zA-Z]+', '', text)

def greek_to_latin(text):
    return ''.join(
        c if 'GREEK' not in unicodedata.name(c) else unicodedata.name(c).split()[-1].title()[:3]
        for c in text
    )

def find_closest_star_name(star_names, text):
    return process.extractOne(text, star_names)[0]

def parse_star_name(input_text, catalogs, star_names):
    input_text = remove_special_characters(input_text)
    input_text = remove_latex_symbols(input_text)
    input_text = greek_to_latin(input_text)

    catalog = find_closest_catalog(catalogs, input_text)
    star_name = find_closest_star_name(star_names, input_text)

    return catalog, star_name

# Example usage:
catalogs = ['HD', 'HR', 'GJ', 'HIP', 'TIC', 'KIC', 'TOI', 'KOI', 'SDSS']
star_names = ['Sirius', 'Wezen', 'Dubhe']

input_text = "H!D_6@7&6\$^!*(Wizier)"  # Example input

catalog, star_name = parse_star_name(input_text, catalogs, star_names)
print(f"Catalog: {catalog}")
print(f"Star Name: {star_name}")
