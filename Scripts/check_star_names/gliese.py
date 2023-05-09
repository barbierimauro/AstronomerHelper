#!/usr/bin/python3
import csv
from fuzzywuzzy import fuzz, process
import statistics


class GlieseSearch:

    def __init__(self, user_string):
        self.user_string = user_string

    def check_number_in_ranges(self, number):
        if ((1 <= float(number) <= 915) or (1001 <= int(number) <= 2057) or (3001 <= int(number) <= 4388) or (9003 <= int(number) <= 9848)):
            return True
        return False

    def normalize_input(self):
        prefixes = ["Gliese", "gliese", "gj", "GJ", "Gl", "GL", "gl", "NN", "Wo"]
        for prefix in prefixes:
            if self.user_string.startswith(prefix):
                self.user_string = self.user_string.replace(prefix, "GJ", 1)
                break
        return self.user_string

    def find_best_matches(self, input_string, names_list):
        methods = [fuzz.partial_ratio, fuzz.token_sort_ratio, fuzz.token_set_ratio, fuzz.ratio]
        results = []

        for name in names_list:
            scores = [method(input_string, name) for method in methods]
            median_score = statistics.median(scores)
            results.append((name, median_score))

        return sorted(results, key=lambda x: x[1], reverse=True)[:10]

    def search(self):
        normalized_string = self.normalize_input()

        if normalized_string:
            number = ''.join(c for c in normalized_string if c.isdigit() or c == '.')
            if not self.check_number_in_ranges(number):
                print("Number is out of range")
                return

            with open("gliese.csv", "r") as csvfile:
                csvreader = csv.DictReader(csvfile)
                correct_names = [row["correct_name"] for row in csvreader]

                if normalized_string in correct_names:
                    print("Exact match found:", normalized_string)
                else:
                    best_matches = self.find_best_matches(normalized_string, correct_names)
                    for match, score in best_matches:
                        print("Match: {} with score: {}".format(match, score))
        else:
            print("Invalid input string")

user_string = input("Enter a string: ").strip()
search = GlieseSearch(user_string)
search.search()
