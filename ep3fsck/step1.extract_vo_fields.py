#!/bin/python3
import csv

# Define the possible keywords
POSSIBLE_KEYWORDS = ['TTYPE', 'TUCD', 'TUNIT', 'TFORM', 'TCOMM', 'TUTYP', 'TDMAX', 'TDMIN', 'TNULL']

# Initialize variables
groups = {}
missing_keywords = {}

# Read the keyword list file
with open('keywords_vo.list', 'r') as f:
    # Loop through the lines of the file
    for line in f:
        # Extract the keyword and group number
        keyword, value = line.strip().split('=')
        keyword = keyword.strip()
        value = value.strip()


        # Check if there is a comment character
        comment_char = '/'
        comment_index = value.find(comment_char)
        if comment_index != -1:
            value = value[:comment_index].strip()

        # Check if the keyword value is enclosed in single quotes
        if value.startswith("'") and value.endswith("'"):
            value = value[1:-1]


        # Determine the group number
        digits = len(keyword) - len(keyword.rstrip('0123456789'))
        group_str = keyword[:len(keyword) - digits]
        group_num = int(keyword[len(keyword) - digits:])
        # Add the keyword to the group
        if group_num not in groups:
            groups[group_num] = {}
        groups[group_num][group_str] = value
        # Check for missing keywords in the group
        if group_num not in missing_keywords:
            missing_keywords[group_num] = set(POSSIBLE_KEYWORDS)
        if group_str in missing_keywords[group_num]:
            missing_keywords[group_num].remove(group_str)

# Count the number of groups
num_groups = len(groups)

# Print missing keywords for each group
for group in groups.keys():
    if group in missing_keywords and missing_keywords[group]:
        missing_keywords_list = ', '.join(missing_keywords[group])
        print(f"COLUMN: {group}, missing keywords: {missing_keywords_list}")
    else:
        print(f"COLUMN: {group}, No missing keywords")

# Create the CSV table
with open('keywords_vo_table.csv', 'w', newline='') as f:
    writer = csv.writer(f)

    # Write the header row
    header = ['COLUMN'] + POSSIBLE_KEYWORDS
    writer.writerow(header)

    # Write the data rows
    for group in groups.keys():
        row = [group]
        for keyword in POSSIBLE_KEYWORDS:
            value = groups[group].get(keyword, '')
            row.append(value)
        writer.writerow(row)

