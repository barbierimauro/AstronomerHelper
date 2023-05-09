#!/usr/bin/python3
# (c) 2023 MAURO BARBIERI   maurobarbieri.science@gmail.com
#
import json
import argparse


def remove_data(node):
    if isinstance(node, dict):
        return {key: remove_data(value) for key, value in node.items()}
    elif isinstance(node, list):
        return [remove_data(value) for value in node]
    else:
        return "__VALUE__"


def print_json_structure(prefix, event, value, level):
    if event == 'map_key':
        print("  " * level + prefix + '.' + value)
        return prefix + '.' + value, level + 1
    elif event == 'start_map':
        return prefix, level
    else:
        return prefix.rpartition('.')[0], level - 1

parser = argparse.ArgumentParser(description="Print the structure of a JSON file using ijson")
parser.add_argument("inputfilename", help="Input JSON file")
args = parser.parse_args()

inputfilename = args.inputfilename

with open(inputfilename, "r") as f:
    data = json.load(f)

structure_data = remove_data(data)
pretty_structure = json.dumps(structure_data, indent=2)

for line in pretty_structure.split('\n'):
    if line.strip() != '"__VALUE__",':
        print(line)
