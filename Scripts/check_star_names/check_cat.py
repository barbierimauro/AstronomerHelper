#!/usr/bin/python3
import re

def is_valid_star_name(star_name, rules):
    for prefix, rule_function in rules.items():
        if star_name.lower().startswith(prefix.lower()):
            if rule_function(star_name):
                return True
    return False

def rewrite_star_name(star_name, rules):
    for prefix, rule_function in rules.items():
        if star_name.lower().startswith(prefix.lower()):
            formatted_name = rule_function(star_name, rewrite=True)
            if formatted_name:
                return formatted_name
    return None

def hd_rules(star_name, rewrite=False):
    pattern = re.compile(r'^(hd|HD|Hd|hD)( {0,1}\d{1,6})([a-zA-Z]{1,2}){0,1}$', re.IGNORECASE)
    match = pattern.match(star_name)
    if match:
        if rewrite:
            formatted_name = match.group(1) + (' ' if not ' ' in match.group(2) else '') + match.group(2).strip() + (match.group(3) if match.group(3) else '')
            return formatted_name.upper()
        else:
            return True
    return False

star_name_rules = {
    "HD": hd_rules,
}

star_names = [
    "HD123",
    "hd 456",
    "Hd789A",
    "hD 321B",
    "HD 12345AB",
    "HD1234",
    "hd98765C",
    "HD 12D",
    "hd3e",
    "HD12345XZ",
    "HR1001",
    "HR 2002",
    "HIP 123",
    "GJ 1001",
]

for star_name in star_names:
    if is_valid_star_name(star_name, star_name_rules):
        print("Valid star name '{}' is rewritten as '{}'".format(star_name, rewrite_star_name(star_name, star_name_rules)))
    else:
        print("Invalid star name '{}'".format(star_name))
