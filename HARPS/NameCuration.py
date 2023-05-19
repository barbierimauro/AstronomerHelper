import pandas as pd
import numpy as np
import re
import sys
import pandas as pd


# Read the CSV file
dataframe = pd.read_csv('harps_all_tables.astro.small.csv', comment='#')

# Group by the 4th column (targets)
grouped = dataframe.groupby(dataframe.columns[3])
# Compute count, mean and standard deviation for 1st and 3rd columns for each group
grouped = grouped.agg(count=('target', 'size'),
                     ra=(dataframe.columns[2], 'mean'),
                     e_ra=(dataframe.columns[2], 'std'),
                     dec=(dataframe.columns[0], 'mean'),
                     e_dec=(dataframe.columns[0], 'std'))

# If standard deviations are NaN or non-defined, replace them with 0.0042
grouped['e_ra'].replace({np.nan: 0.0042}, inplace=True)
grouped['e_dec'].replace({np.nan: 0.0042}, inplace=True)

# Reset index
grouped.reset_index(inplace=True)


# Sort DataFrame by 'Count'
grouped.sort_values(by='count', ascending=False, inplace=True)

# Write to CSV
grouped.to_csv('targets.csv', index=False)

constellations = ["AND", "ANT", "APS", "AQR", "AQL", "ARA", "ARI", "AUR", "BOO", "CAE", "CAM", "CNC", "CVN", "CMA", "CMI", "CAP", "CAR", "CAS", "CEN", "CEP", "CET", "CHA", "CIR", "COL", "COM", "CRA", "CRB", "CRV", "CRT", "CRU", "CYG", "DEL", "DOR", "DRA", "EQU", "ERI", "FOR", "GEM", "GRU", "HER", "HOR", "HYA", "HYI", "IND", "LAC", "LEO", "LMI", "LEP", "LIB", "LUP", "LYN", "LYR", "MEN", "MIC", "MON", "MUS", "NOR", "OCT", "OPH", "ORI", "PAV", "PEG", "PER", "PHE", "PIC", "PSC", "PSA", "PUP", "PYX", "RET", "SGE", "SGR", "SCO", "SCL", "SCT", "SER", "SEX", "TAU", "TEL", "TRI", "TRA", "TUC", "UMA", "UMI", "VEL", "VIR", "VOL", "VUL"]
greekletters = ["ALF", "BET", "GAM", "DEL", "EPS", "ZET", "ETA", "THE", "IOT", "KAP", "LAM", "MU", "NU", "XI", "OMIC", "PI", "RHO", "SIG", "TAU", "UPS", "PHI", "CHI", "PSI", "OME"]
targets_to_modify = ["0135_BB_SCL","0303_EW_ERI","0314_TZ-FOR","0411_GW-ERI","0446_AL-DOR","0530_TZ-MEN","0530_UX-MEN","0628_V722-MEN","0630_KL-CMA","0654_AUMON","0658_V358-PUP","0749_V397-PUP","0826_NO-PUP","0831_VZ-HYA","0910_PT-VEL","0921_NY-HYA","1934_V4089_SGR","2028_MP_DEL","2231_KX-AQR","2234_LL-AQR"]
solsys_list = ['ASTRAEA', 'CERES', 'COMET', 'EUROPA', 'GANYMEDE', 'GANIMEDE', 'GIOVE', 'JUNO', 'JUPITER', 'IRENE', 'IRIS', 'HERCULINA', 'MASSALIA', 'MELPOMENE', 'PALLAS', 'PARTHENOPE', 'TITAN', 'VENUS', 'VESTA', 'MOON', 'SKY']
replacement_dict = {'6-HEBE':'HEBE', 'MOOON':'MOON', 'ZENITH':'SKY', 'SOLAR':'SUN', 'SUN':'SUN', 'TWILIGHT':'SKY'}


# Define a function to perform the specific operation
def curate_target(target):
    target_curated=target
    target_is_curated = False
    curation = -1

    if target == "18_SCO":
        print(target)

    if not target_is_curated and any(solsys in target for solsys in solsys_list):
        target_curated = [solsys for solsys in solsys_list if solsys in target][0]
        curation = 100001
        target_is_curated = True
    elif not target_is_curated and target == 'IO':
        target_curated = 'IO'
        curation = 100001
        target_is_curated = True
    elif not target_is_curated and any(solsys in target for solsys in replacement_dict.keys()) :
        target_curated = replacement_dict[[solsys for solsys in replacement_dict.keys() if solsys in target][0]]
        curation = 100001
        target_is_curated = True
    elif target.startswith('HD')  and not target_is_curated:
        target_curated = target_curated.replace('-', '').replace(' ', '')
        parts = target_curated.split('_')
        if len(parts) > 2:
            target_curated = parts[0] + parts[1]
        target_curated = target_curated.replace('_', '')
        curation = 1
        target_is_curated = True
    elif target.startswith(('GJ', 'GL', 'HIP', 'HR')) and not target_is_curated:
        target_curated = target_curated.replace('-', '').replace('_', '').replace(' ', '')
        curation = 2
        target_is_curated = True
    elif target.startswith('TOI') and not target_is_curated:
        target_curated = target_curated.replace('-', ' ').replace('_', ' ')
        curation = 3
        target_is_curated = True
    elif target.startswith('EPIC') and not target_is_curated:
        target_curated = target_curated.replace('-', ' ')
        curation = 4
        target_is_curated = True
    elif target.startswith('WASP') and not target_is_curated:
        target_curated = target_curated.replace('_ORB', '').replace('_TRA', '')
        curation = 5
        target_is_curated = True
    elif target.startswith('TYC') and not target_is_curated:
        target_curated = target_curated.replace('_', ' ').replace('TYC-', 'TYC ')
        curation = 6
        target_is_curated = True
    elif target.startswith('BD') and not target_is_curated:
        target_curated = target_curated.replace('_', ' ').replace('--', '-').replace('BDCAP', 'BD CAP')
        curation = 7
        target_is_curated = True
    elif target.startswith('$EPS') and not target_is_curated:
        target_curated = target_curated.replace('$EPSILON$,CMA', 'EPS CMA')
        curation = 8
        target_is_curated = True
    elif target.startswith('$IOT') and not target_is_curated:
        target_curated = target_curated.replace('$IOTA$,CAR', 'IOT CAR')
        curation = 9
        target_is_curated = True
    elif target.startswith('ALF') and not target_is_curated:
        target_curated = target_curated.replace('-', ' ').replace('_', ' ')
        curation = 10
        target_is_curated = True
    elif target.startswith('ALP') and not target_is_curated:
        target_curated = target_curated.replace('-', ' ').replace('_', ' ').replace('CENA', 'CEN A').replace('CENB', 'CEN B').replace('CENTAURI', 'CEN').replace('ALPHA', 'ALF ').replace('ALP', 'ALF ')
        curation = 11
        target_is_curated = True
    elif target.startswith('BET') and not target_is_curated:
        target_curated = target_curated.replace('BETA', 'BET ').replace('-', ' ').replace('_', ' ').replace('BETCEN', 'BET CEN').replace('BETCRT', 'BET CRT').replace('BETHYI', 'BET HYI')
        curation = 12
        target_is_curated = True
    elif target.startswith('GAM') and not target_is_curated:
        target_curated = target_curated.replace('GAMMA', 'GAM ').replace('-', ' ').replace('_', ' ').replace('GAMCIR', 'GAM CIR').replace('GAMCEN', 'GAM CEN')
        curation = 13
        target_is_curated = True
    elif target.startswith('DEL') and not target_is_curated:
        target_curated = target_curated.replace('DELTA', 'DEL ').replace('-', ' ').replace('_', ' ').replace('DELERI', 'DEL ERI')
        curation = 14
        target_is_curated = True
    elif target.startswith('EPS') and not target_is_curated:
        target_curated = target_curated.replace('EPSILON', 'EPS ').replace('-', ' ').replace('_', ' ').replace('EPSERI', 'EPS ERI').replace('EPSLUP', 'EPS LUP').replace('EPSIND', 'EPS IND').replace('ERIDANI', 'ERI').replace('EPS-INDI-A', 'EPS IND')
        curation = 15
        target_is_curated = True
    elif target.startswith('TAU') and not target_is_curated:
        target_curated = target_curated.replace('-', ' ').replace('_', ' ').replace('CETI', 'CET').replace('TAUCET', 'TAU CET').replace('TAUBOO', 'TAU BOO')
        curation = 16
        target_is_curated = True
    elif target.startswith('IOT') and not target_is_curated:
        target_curated = target_curated.replace('-', ' ').replace('_', ' ').replace('IOTHOR', 'IOT HOR').replace('IOTCAR', 'IOT CAR')
        curation = 17
        target_is_curated = True
    elif target.startswith('ZET') and not target_is_curated:
        target_curated = target_curated.replace('ZETA', 'ZET ').replace('-', ' ').replace('_', ' ').replace('ZET02RET', 'ZET02 RET').replace('ZET1RET', 'ZET1 RET').replace('ZET2RET', 'ZET2 RET').replace('ZETTUC', 'ZET TUC')
        curation = 18
        target_is_curated = True
    elif target.startswith(('LRC0', 'LRA0', 'IRA0', 'SRC0', 'SRA0')) and not target_is_curated:
        target_curated = target_curated.replace('_', ' ')
        curation = 19
        target_is_curated = True
    elif target.startswith('LHS') and not target_is_curated:
        target_curated = target_curated.replace('LHS', 'LHS ')
        if '_' in target_curated:
            target_curated = target_curated[:target_curated.index('_')]
            curation = 20
            target_is_curated = True
    elif target.startswith('CD') and not target_is_curated:
        target_curated = target_curated.replace('--', '-').replace('CD-53-251', 'CD-53 251').replace('CD-69 1055', 'CD-69 1055').replace('CD-71-1234', 'CD-71 1234').replace('CD-84-0080', 'CD-84 0080').replace('CD-TAU', 'CD TAU').replace('CD_CIR', 'CD CIR')
        curation = 21
        target_is_curated = True
    elif target.startswith('HV') and not target_is_curated:
        target_curated = target_curated.replace('HV-LUP', 'HV LUP').replace('HV_LUP', 'HV LUP').replace('HVERI', 'HV ERI')
        curation = 22
        target_is_curated = True
    elif target.startswith('PROXIMA') and not target_is_curated:
        target_curated = 'HIP 70890'
        curation = 23
        target_is_curated = True
    elif target.startswith('ACHERNAR') and not target_is_curated:
        target_curated = 'HD 10144'
        curation = 24
        target_is_curated = True
    elif target.startswith('PROCYON') and not target_is_curated:
        target_curated = 'HD 61421'
        curation = 25
        target_is_curated = True
    elif target.startswith('LP') and not target_is_curated:
        target_curated = target_curated.replace('LP771-95A', 'BD-17 588A').replace('LP_LIB', 'LP LIB').replace('LP_MUS', 'LP MUS').replace('LP_VIR', 'LP VIR')
        curation = 26
        target_is_curated = True
    elif target.startswith('NU') and not target_is_curated:
        target_curated = target_curated.replace('NU-IND', 'NU IND').replace('NUPHE', 'NU PHE').replace('NU_OCT', 'NU OCT')
        curation = 27
        target_is_curated = True
    elif target.startswith('OMI') and not target_is_curated:
        target_curated = target_curated.replace('-', ' ').replace('OMILUP', 'OMI LUP')
        curation = 28
        target_is_curated = True
    elif target.startswith('NLTT') and not target_is_curated:
        target_curated = target_curated.replace('-', ' ').replace('NLTT44619_F2V', 'NLTT44619')
        curation = 29
        target_is_curated = True
    elif target.startswith('CPD') and not target_is_curated:
        target_curated = target_curated.replace('CPD+60-944A', 'CPD-60 944A').replace('CPD-28-2561', 'CPD-28 2561').replace('CPD-43-7188', 'CPD-43 7188').replace('CPD-60-944A', 'CPD-60 944A')
        curation = 30
        target_is_curated = True
    elif target in targets_to_modify and not target_is_curated:
        if target == "0654_AUMON":
            target_curated = "AU MON"
            curation = 31
            target_is_curated = True
        elif target == "0108 PHE":
            target_curated = "108 PHE"
            curation = 32
            target_is_curated = True
        else:
            target_curated = target[5:].replace("_", " ").replace("-", " ")
            curation = 33
            target_is_curated = True
    elif any(target.endswith(c) for c in constellations) and not target_is_curated:
        if target.startswith("V ") or target.startswith("V*"):
            target_curated = target[2:]
            curation = 42
            target_is_curated = True
        elif target[0] == "V" and target[1].isdigit():
            constellation = target[-3:]
            remaining = target[:-3].replace("-", " ").replace("_", " ")
            target_curated = remaining + " " + constellation
            curation = 43
            target_is_curated = True
    elif  target_is_curated == False and any(target.endswith(c) for c in constellations):
        for constellation in constellations:
            if target.endswith(constellation):
                prefix = target[:-len(constellation)].strip()  # the part before constellation
                if len(prefix) == 2 and prefix.isalpha() and not target_is_curated:  
                    target_curated = prefix + " " + constellation
                    curation = 34
                    target_is_curated = True
                elif len(prefix) == 1 and prefix.isalpha() and not target_is_curated:  
                    target_curated = prefix + " " + constellation
                    curation = 35
                    target_is_curated = True
                elif len(prefix) > 1 and prefix[0].isalpha() and not target_is_curated:  
                    target_curated = prefix[0] + " " + constellation
                    curation = 36
                    target_is_curated = True
                elif len(prefix) > 2 and prefix[:2].isalpha() and not target_is_curated:
                    if any(char.isdigit() or not char.isalpha() for char in prefix[2:]):
                        target_curated = prefix[:2] + " " + constellation
                        curation = 37
                        target_is_curated = True
                else:
                    for i in range(3, 0, -1):  # this loop will check for 3, 2, and 1 digit(s) at the start
                        if len(prefix) >= i and prefix[:i].isdigit() and not target_is_curated:
                            target_curated = prefix[:i] + " " + constellation
                            curation = 38
                            target_is_curated = True
                            break
    elif target.startswith('NGC') and not target_is_curated:
        if target.startswith('NGC '):
            target_curated = target
            curation = 39
            target_is_curated = True
        elif target[3].isdigit():
            target_curated = "NGC " + target[3:7] + " "
            target_curated += ''.join(char for char in target[7:] if char.isdigit())
            curation = 40
            target_is_curated = True
    elif target.startswith('IC') and not target_is_curated:
        if target[2].isdigit():
            target_curated = "IC " + target[2:6] + " "
            target_curated += ''.join(char for char in target[6:] if char.isdigit())
            curation = 41
            target_is_curated = True
    else:
        target_curated = "NONCURATED"
        curation = -1
        target_is_curated = True

    return target_curated, curation


# create an empty list
curated_targets_list = []

# iterate over the rows in the grouped DataFrame
for row in grouped.itertuples():
    # curate the target and get curation status
    curated_target, curation = curate_target(row.target)

    # Create a new row with all original columns, curated target, and curation status
    new_row = row._asdict()  # Convert row to OrderedDict
    new_row['target_curated'] = curated_target
    new_row['curation'] = curation
    # Add the new row to the list
    curated_targets_list.append(new_row)

# Create a new DataFrame from the list
curated_targets = pd.DataFrame(curated_targets_list)



# Write to CSV
curated_targets.to_csv('target_curated.csv', index=False)


#print(curated_targets)


sys.exit()

cata_list = ["HD", "HR", "HIP", "BD", "CD", "CPD", "GJ", "GL", "TIC", "WASP", "TYC", "2MASS", "K2", "COROT", "EPIC", "HAT", "HV", "IC", "MELOTTE", "KELT", "LHS", "LMC", "LP", "IRA", "LRA", "LRC", "SRA", "SRC", "NGC", "NLTT", "LTT", "ROSS", "NOI", "OGLE", "QATAR", "SERAM", "UCAC2", "UCAC4", "SAND", "SAO", "RX", "SMC", "SW", "TOI", "SUN", "MOON", "ALP", "ALF", "BET", "GAM", "DEL", "EPS", "ZET", "ETA", "IOT", "KAP", "LAM", "MU", "NU", "OMI", "PI", "RHO", "SIG", "TAU", "UPS", "PHI", "CHI", "PSI", "OME", "PROXIMA", "SIRIUS", "PROCYON", "ACHERNAR",
"SOLAR", "SKY", "TEST", "B_", "\$EPS" , "\$IOT", "GANI", "GANY", "EUROPA", "VESTA", "CERES", "VENUS", "JUPITER", "GIOVE", "NAME PROXIMA", "NAME BARNARD"]

dataframe = pd.read_csv('harps_all_tables.astro.small.csv', comment='#')

# Preprocess the dataframe by replacing "$," with "$ " and removing all "$"
dataframe = dataframe.replace({"$,": "$ ", "$": ""}, regex=True)

total_count = 0
counts = {}
for cata in cata_list:
    count = dataframe.iloc[:, 3].str.contains("^" + cata, case=False, regex=True).sum()
    total_count += count
    counts[cata] = count
    #print(f"{cata} {count}")

# Sort dictionary by value
sorted_counts = dict(sorted(counts.items(), key=lambda item: item[1]))

for cata, count in sorted_counts.items():
    print(f"{cata} {count}")

print(f"total {total_count}")

# Removing all rows where the 3rd column's values start with any string from the cata_list
#pattern = '^' + '|^'.join(cata_list)
#dataframe = dataframe[~dataframe.iloc[:, 3].str.contains(pattern, case=False, regex=True)]

#dataframe.to_csv('remnant.csv', index=False)



sys.exit()


sys.exit()

# Read the CSV file
dataframe = pd.read_csv('harps_all_tables.astro.small.csv', comment='#')

# Group by the 4th column (targets)
grouped = dataframe.groupby(dataframe.columns[3])

# Compute count, mean and standard deviation for 1st and 3rd columns for each group
result = grouped.agg(count=('target', 'size'),
                     ra=(dataframe.columns[2], 'mean'),
                     e_ra=(dataframe.columns[2], 'std'),
                     dec=(dataframe.columns[0], 'mean'),
                     e_dec=(dataframe.columns[0], 'std'))

# If standard deviations are NaN or non-defined, replace them with 0.0042
result['e_ra'].replace({np.nan: 0.0042}, inplace=True)
result['e_dec'].replace({np.nan: 0.0042}, inplace=True)

# Reset index
result.reset_index(inplace=True)


# Sort DataFrame by 'Count'
result.sort_values(by='count', ascending=False, inplace=True)

# Write to CSV
result.to_csv('targets.csv', index=False)


sys.exit()
