#!/usr/bin/python3
# (c) 2023 MAURO BARBIERI   maurobarbieri.science@gmail.com
#
# clean the JSON shit and provide a human readable output
#
"""
This script reads a JSON file and prints the keywords (keys) present in the file along with the number of times they occur.

Usage:
python3 keywords_counter.py <inputfilename>

Arguments:
inputfilename: Path to the input JSON file.

Example:
python3 keywords_counter.py example.json

Functions:
count_keywords(data, counter):
Recursively traverses the given JSON data structure and updates the keyword occurrences count.

Args:
    data (Union[dict, list]): JSON data structure.
    counter (collections.Counter): A counter object to keep track of the keyword occurrences.

Returns:
    None
"""

import json
import argparse
from collections import Counter

# Define a function to recursively count the occurrences of keywords in the JSON data
def count_keywords(data, counter):
    # If the data is a dictionary, iterate through its key-value pairs
    if isinstance(data, dict):
        for key, value in data.items():
            # Increment the count for the current key
            counter[key] += 1
            # Recursively process the value
            count_keywords(value, counter)
    # If the data is a list, iterate through its items
    elif isinstance(data, list):
        for item in data:
            # Recursively process the item
            count_keywords(item, counter)
    return

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Read a JSON file and print the keywords and their occurrences")
parser.add_argument("inputfilename", help="Input JSON file")
args = parser.parse_args()
inputfilename = args.inputfilename

# Read the JSON file
with open(inputfilename) as f:
    data = json.load(f)

# Create a counter to keep track of keyword occurrences
counter = Counter()
# Count the keywords in the JSON data
count_keywords(data, counter)
# Print the sorted keywords and their counts
for keyword, count in sorted(counter.items(), key=lambda x: x[1], reverse=True):
    print(f"{keyword}: {count}")
