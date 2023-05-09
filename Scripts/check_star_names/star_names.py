#!/usr/bin/python3
import sys
import os
import csv
import pandas as pd
from fuzzywuzzy import fuzz
from numpy import median


constellations = [
    'Andromeda', 'Antlia', 'Apus', 'Aquarius', 'Aquila', 'Ara', 'Aries', 'Auriga', 'Bootes',
    'Caelum', 'Camelopardalis', 'Cancer', 'Canes Venatici', 'Canis Major', 'Canis Minor',
    'Capricornus', 'Carina', 'Cassiopeia', 'Centaurus', 'Cepheus', 'Cetus', 'Chamaeleon',
    'Circinus', 'Columba', 'Coma Berenices', 'Corona Australis', 'Corona Borealis', 'Corvus',
    'Crater', 'Crux', 'Cygnus', 'Delphinus', 'Dorado', 'Draco', 'Equuleus', 'Eridanus',
    'Fornax', 'Gemini', 'Grus', 'Hercules', 'Horologium', 'Hydra', 'Hydrus', 'Indus', 'Lacerta',
    'Leo', 'Leo Minor', 'Lepus', 'Libra', 'Lupus', 'Lynx', 'Lyra', 'Mensa', 'Microscopium',
    'Monoceros', 'Musca', 'Norma', 'Octans', 'Ophiuchus', 'Orion', 'Pavo', 'Pegasus', 'Perseus',
    'Phoenix', 'Pictor', 'Pisces', 'Piscis Austrinus', 'Puppis', 'Pyxis', 'Reticulum', 'Sagitta',
    'Sagittarius', 'Scorpius', 'Sculptor', 'Scutum', 'Serpens', 'Sextans', 'Taurus', 'Telescopium',
    'Triangulum', 'Triangulum Australe', 'Tucana', 'Ursa Major', 'Ursa Minor', 'Vela', 'Virgo',
    'Volans', 'Vulpecula'
]

constellations_abbreviations = [
    'And', 'Ant', 'Aps', 'Aqr', 'Aql', 'Ara', 'Ari', 'Aur', 'Boo',
    'Cae', 'Cam', 'Cnc', 'CVn', 'CMa', 'CMi',
    'Cap', 'Car', 'Cas', 'Cen', 'Cep', 'Cet', 'Cha',
    'Cir', 'Col', 'Com', 'CrA', 'CrB', 'Crv',
    'Crt', 'Cru', 'Cyg', 'Del', 'Dor', 'Dra', 'Equ', 'Eri',
    'For', 'Gem', 'Gru', 'Her', 'Hor', 'Hya', 'Hyi', 'Ind', 'Lac',
    'Leo', 'LMi', 'Lep', 'Lib', 'Lup', 'Lyn', 'Lyr', 'Men', 'Mic',
    'Mon', 'Mus', 'Nor', 'Oct', 'Oph', 'Ori', 'Pav', 'Peg', 'Per',
    'Phe', 'Pic', 'Psc', 'PsA', 'Pup', 'Pyx', 'Ret', 'Sge',
    'Sgr', 'Sco', 'Scl', 'Sct', 'Ser', 'Sex', 'Tau', 'Tel',
    'Tri', 'TrA', 'Tuc', 'UMa', 'UMi', 'Vel', 'Vir',
    'Vol', 'Vul'
]

constellations_genitive = [
    'Andromedae', 'Antliae', 'Apodis', 'Aquarii', 'Aquilae', 'Arae', 'Arietis', 'Aurigae', 'Boötis',
    'Caeli', 'Camelopardalis', 'Cancri', 'Canum Venaticorum', 'Canis Majoris', 'Canis Minoris',
    'Capricorni', 'Carinae', 'Cassiopeiae', 'Centauri', 'Cephei', 'Ceti', 'Chamaeleontis',
    'Circini', 'Columbae', 'Comae Berenices', 'Coronae Australis', 'Coronae Borealis', 'Corvi',
    'Crateris', 'Crucis', 'Cygni', 'Delphini', 'Doradus', 'Draconis', 'Equulei', 'Eridani',
    'Fornacis', 'Geminorum', 'Gruis', 'Herculis', 'Horologii', 'Hydrae', 'Hydri', 'Indi', 'Lacertae',
    'Leonis', 'Leonis Minoris', 'Leporis', 'Librae', 'Lupi', 'Lyncis', 'Lyrae', 'Mensae', 'Microscopii',
    'Monocerotis', 'Muscae', 'Normae', 'Octantis', 'Ophiuchi', 'Orionis', 'Pavonis', 'Pegasi', 'Persei',
    'Phoenicis', 'Pictoris', 'Piscium', 'Piscis Austrini', 'Puppis', 'Pyxidis', 'Reticuli', 'Sagittae',
    'Sagittarii', 'Scorpii', 'Sculptoris', 'Scuti', 'Serpentis', 'Sextantis', 'Tauri', 'Telescopii',
    'Trianguli', 'Trianguli Australis', 'Tucanae', 'Ursae Majoris', 'Ursae Minoris', 'Velorum', 'Virginis',
    'Volantis', 'Vulpeculae'
]

greek_letters = [
    'Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Zeta', 'Eta', 'Theta', 'Iota', 'Kappa', 'Lambda',
    'Mu', 'Nu', 'Xi', 'Omicron', 'Pi', 'Rho', 'Sigma', 'Tau', 'Upsilon', 'Phi', 'Chi', 'Psi', 'Omega'
]

# Convert each letter to its abbreviation
greek_letters_abbreviations = [x[:3] for x in greek_letters]


CACHE_DIR = os.path.expanduser('~/.cache')
os.makedirs(CACHE_DIR, exist_ok=True)
IAUCSN_FILENAME = "IAU-CSN.csv"
iaucsn_file_path = os.path.join(CACHE_DIR, IAUCSN_FILENAME)

if not os.path.exists(iaucsn_file_path):
    url = "https://raw.githubusercontent.com/barbierimauro/ESOarchive/main/data/IAU-CSN.csv"
    iaucsn = pd.read_csv(url)
    iaucsn.to_csv(iaucsn_file_path, index=False)
else:
    iaucsn = pd.read_csv(iaucsn_file_path)


# Define constants
score_cutoff = 60
LENGTH_PENALTY_FACTOR = 10
SPACE_PENALTY = 100

# 1) Record scores of each fuzzy matching method separately
def calculate_scores(input_str, candidate_str):
    methods = [fuzz.partial_ratio, fuzz.token_sort_ratio, fuzz.token_set_ratio, fuzz.ratio]
    scores = [method(input_str.lower(), candidate_str.lower()) for method in methods]
    return scores


# 2) Exact matches have scores of 1000
def exact_match_bonus(input_str, candidate_str, scores):
    if input_str.lower() == candidate_str.lower():
        return [1000] * len(scores)
    return scores


# 3) Add bonus scores for character sequence matches
def sequence_bonus(input_str, candidate_str, scores):
    for n in range(6, 3, -1):
        if any(input_str.lower()[i:i+n] in candidate_str.lower() for i in range(len(input_str)-n+1)):
            bonus = {6: 50, 5: 20, 4: 5}[n]
            return [score + bonus for score in scores]
    return scores


# 4) Calculate the median of the scores
def median_score(scores):
    return median(scores)


# Apply penalties
def apply_penalties(input_str, candidate_str, score):
    length_ratio = abs(len(input_str) - len(candidate_str)) / max(len(input_str), len(candidate_str))
    length_penalty = LENGTH_PENALTY_FACTOR * length_ratio
    space_penalty = SPACE_PENALTY if ' ' in candidate_str else 0
    final_score = score - length_penalty - space_penalty
    return final_score




# Initialize data
STARINI = input("Enter object name: ")
STARINI = ' '.join(STARINI.split())  # Remove multiple spaces and standardize the input

results = []

# Perform fuzzy search of STARINI using names in iaucsn dataframe
for index, row in iaucsn.iterrows():

    scores = calculate_scores(STARINI, row['Name'])
    scores = exact_match_bonus(STARINI, row['Name'], scores)
    scores = sequence_bonus(STARINI, row['Name'], scores)
    median_score_value = median_score(scores)

    final_score = apply_penalties(STARINI.lower(), row['Name'].lower(), median_score_value)
    if final_score >= score_cutoff:
        results.append({"Name": row["Name"], "median_score": final_score,  "Scores": scores})




# Check against combinations of Greek letters and constellations
for greek_letter in greek_letters + greek_letters_abbreviations:
    for constellation in constellations + constellations_genitive + constellations_abbreviations:
        combined_name = "{} {}".format(greek_letter, constellation)

        scores = calculate_scores(STARINI, combined_name)
        scores = exact_match_bonus(STARINI, combined_name, scores)
        scores = sequence_bonus(STARINI, combined_name, scores)

        median_score_value = median_score(scores)

        final_score = apply_penalties(STARINI.lower(), combined_name.lower(), median_score_value)
        #print(scores)
        if final_score >= score_cutoff:
            results.append({"Name": combined_name, "median_score": final_score, "Scores": scores})



# Sort results by the median of their scores
sorted_results = sorted(results, key=lambda x: x['median_score'], reverse=True)

# Print results
if len(sorted_results) == 0:
    print("No results found.")
else:
    print("{} results found. Displaying top 20:".format(len(sorted_results)))
    for i, result in enumerate(sorted_results[:20]):
        print(i,result)








