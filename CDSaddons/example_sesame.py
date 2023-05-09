#!/bin/python3
from sesame_query import SesameQuery

# Instantiate the SesameQuery class
myquery = SesameQuery()

# Search for the object "sigma oct"
result = myquery.search("HD-123")

# Print the complete results
for key, value in result.items():
    print(f"{key}: {value}")
