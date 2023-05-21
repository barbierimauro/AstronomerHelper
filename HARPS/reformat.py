import pandas as pd
import csv
import re

# Load the data
df = pd.read_csv('targets.csv')

# Specific replacements
replacements = {
    'TCHOUBIDOO1':'WASP-17',
    'TCHOUBIDOO2':'WASP-15',
    'TCHOUBIDOO3':'WASP-15',
    "$EPSILON$,CMA": 'EPS CMA',
    '$IOTA$,CAR': 'IOT CAR',
    'HD,151804': 'HD 151804',
    'I,CAR': 'V* I CAR',
    'BD+103022P':'TYC 964-1658-1',
    'BD+103022S':'TYC 964-160-1',
    'BD-013612P':'CSI-01 3612 1',
    'BD-013612S':'CSI-01 3612 2',
    'B-CD71-1579': 'CD-71 1579',
    'B-CD-80-28': 'CD-80 28',
    'B-TYC-8620': 'TYC 8620',
    'B-V936-CEN': 'V* V936 CEN',
    'AHD187669': 'HD187669',
    'AHD231630': 'HD231630',
    'BHD104856': 'HD104856',
    'BHD312256': 'HD312256',
    'BHD327692': 'HD327692',
    'BDCAP': 'V* BD CAP',
    'BD_12_4982_9.3' : 'BD-12 4982',
    'BD20594_PLEI':'BD+20 594',
    'LP771-95A':'BD-17 588A',
    'LP_LIB': 'V* LP LIB',
    'LP_MUS': 'V* LP MUS',
    'LP_VIR': 'V* LP VIR',
    'LPV-41682':'OGLE LMC-LPV-41682',
    'LPV41682':'OGLE LMC-LPV-41682',
    'TIC92304420 - 03H22':'TIC92304420',
    'HR4963SPHOT':'HR 4963',
    'HR CAR':'V* HR CAR',
    'HRCAR':'V* HR CAR',
    'HR-CAR':'V* HR CAR',
    'GL-CAR':'V* GL CAR',
    'GI581' : 'GJ581',
    'ROSS 128':'ROSS 128',
    'ACHERNAR':'ACHERNAR',
    'PROCYON':'PROCYON',
    'SIRIUS':'SIRIUS',
    'CANOPUS':'CANOPUS',
    'ACRUX':'ACRUX',
    'ARCTURUS':'ARCTURUS',
    'SIRIUS':'SIRIUS',
    'FEIGE86':'FEIGE86',
    'SW1143+0633UPPER':'SW1143+0633',
    'TWHYDRA':'TW HYA',
    'CCDM-J10573-6902A-B':'CCDM J10573-6902AB',
    'CCDM-J12377-2708AB_F1V':'CCDM J12377-2708AB',
    'PEACOCK':'PEACOCK',
    '0654_AUMON':'V* AU MON',
    '0108_PHE':'ZET PHE',
    '0400_TAU':'LAM TAU',
    'A-CAR':'A CAR',
    'B-CAR':'B CAR',
    'C-SCO':'C SCO',
    'A-VEL':'A VEL',
    'B-VEL':'B VEL',
    'C-VEL':'C VEL',
#    '':'',
}
# notes to the list
# "HR CAR" for protect it and do not modify with other rules
# "GL CAR" for protect it and do not modify with other rules
# add also 
# ACHERNAR PROXIMA PROCYON BARNARD FEIGE86 SIRIUS
# CANOPUS ACRUX ARCTURUS

GLIESE_list = ['GJ','GJ','GLIESE']
constellations = ["AND", "ANT", "APS", "AQR", "AQL", "ARA", "ARI", "AUR", "BOO", "CAE", 
                  "CAM", "CNC", "CVN", "CMA", "CMI", "CAP", "CAR", "CAS", "CEN", "CEP", 
                  "CET", "CHA", "CIR", "COL", "COM", "CRA", "CRB", "CRV", "CRT", "CRU", 
                  "CYG", "DEL", "DOR", "DRA", "EQU", "ERI", "FOR", "GEM", "GRU", "HER", 
                  "HOR", "HYA", "HYI", "IND", "LAC", "LEO", "LMI", "LEP", "LIB", "LUP", 
                  "LYN", "LYR", "MEN", "MIC", "MON", "MUS", "NOR", "OCT", "OPH", "ORI", 
                  "PAV", "PEG", "PER", "PHE", "PIC", "PSC", "PSA", "PUP", "PYX", "RET", 
                  "SGE", "SGR", "SCO", "SCL", "SCT", "SER", "SEX", "TAU", "TEL", "TRI", 
                  "TRA", "TUC", "UMA", "UMI", "VEL", "VIR", "VOL", "VUL"]
greekletters = ["ALF", "BET", "GAM", "DEL", "EPS", "ZET", "ETA", "THE", "IOT", "KAP", 
                "LAM", "MU" , "NU" , "XI" , "OMI", "PI" , "RHO", "SIG", "TAU", "UPS", 
                "PHI", "CHI", "PSI", "OME"]
solsys_list = ['ASTRAEA', 'CERES', 'COMET', 'EUROPA', 'GANYMEDE', 'JUNO', 'JUPITER', 'IRENE', 'IRIS', 'HERCULINA', 
               'MASSALIA', 'MELPOMENE', 'PALLAS', 'PARTHENOPE', 'TITAN', 'VENUS', 'VESTA', 'MOON', 'SKY']
solsys_repl = {'GIOVE':'JUPITER','GANIMEDE':'GANYMEDE', '6-HEBE':'HEBE', 'MOOON':'MOON','THE-MOON':'MOON', 'ZENITH':'SKY', 'SOLAR':'SKY', 'SUN':'SKY', 'TWILIGHT':'SKY'}


# Define a function to perform the transformations
def curate_name(name):
    original_name = name
    curation=0  # real part
    kuration=0  # imaginary part
    note = ""
    target_is_curated = False
    ncur = 0

    # pre-processing
    name=name.replace('--', ' ')
    # Collapse multiple spaces into one
    name = ' '.join(name.split())
    # basic curation
    if name[0] in ['A', 'B', 'C', 'D', 'E', 'F'] and name[0].isalpha() and name[1] in ['_', '-'] and name[2:4] in ['HD', 'HR']:
        name = name[2:]
        kuration = -1
    if name[0] in ['A', 'B', 'C', 'D', 'E', 'F'] and name[0].isalpha() and name[1] in ['_', '-'] and name[2:5] in ['HIP', 'TYC']:
        name = name[2:]
        kuration = -10+kuration
    if name.startswith(('HD_', 'HD-')):
        name = name.replace('_', '',1).replace('-', '',1)
        kuration = -100+kuration
    if name.startswith(('HR_', 'HR-')):
        name = name.replace('_', '',1).replace('-', '',1)
        kuration = -1000+kuration
    if name.startswith(('HIP_', 'HIP-')):
        name = name.replace('_', '',1).replace('-', '',1)
        kuration = -10000+kuration
    if name.startswith('GC-'):
        name = name.replace('GC-', '')
        kuration = -100000+kuration

    #solar system
    if original_name == 'IO' and not target_is_curated:
        name = 'IO'
        curation = 100001
        target_is_curated = True
        ncur = ncur +1
    if any(solsys in name for solsys in solsys_list) and not target_is_curated:
        name = [solsys for solsys in solsys_list if solsys in name][0]
        curation = 100001
        target_is_curated = True
        ncur = ncur +1
    if any(solsys in name for solsys in solsys_repl.keys()) and not target_is_curated:
        name = solsys_repl[[solsys for solsys in solsys_repl.keys() if solsys in name][0]]
        curation = 100001
        target_is_curated = True
        ncur = ncur +1

    # Stars
    if name in replacements and not target_is_curated:
        name = replacements[name]
        target_is_curated = True
        curation = 1
        ncur = ncur +1
    if 'PROXIMA' in original_name and not target_is_curated:
        name = 'HIP 70890'
        target_is_curated = True
        curation = 1
        ncur = ncur +1
    if 'BARNARD' in original_name and not target_is_curated:
        name = 'HIP 87937'
        target_is_curated = True
        curation = 1
        ncur = ncur +1
    # BAYER
    if name.startswith('ALF') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ')
        curation = 101
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('ALP') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ').replace('CENA', 'CEN A').replace('CENB', 'CEN B').replace('CENTAURI', 'CEN').replace('ALPHA', 'ALF ').replace('ALP', 'ALF ')
        curation = 102
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('BET') and not target_is_curated:
        name = name.replace('BETA', 'BET ').replace('-', ' ').replace('_', ' ').replace('BETCEN', 'BET CEN').replace('BETCRT', 'BET CRT').replace('BETHYI', 'BET HYI')
        curation = 103
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('GAM') and not target_is_curated:
        name = name.replace('GAMMA', 'GAM ').replace('-', ' ').replace('_', ' ').replace('GAMCIR', 'GAM CIR').replace('GAMCEN', 'GAM CEN')
        curation = 104
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('DEL') and not target_is_curated:
        name = name.replace('DELTA', 'DEL ').replace('-', ' ').replace('_', ' ').replace('DELERI', 'DEL ERI')
        curation = 105
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('EPS') and not target_is_curated:
        name = name.replace('EPSILON', 'EPS ').replace('-', ' ').replace('_', ' ').replace('EPSERI', 'EPS ERI').replace('EPSLUP', 'EPS LUP').replace('EPSIND', 'EPS IND').replace('ERIDANI', 'ERI').replace('EPS-INDI-A', 'EPS IND')
        curation = 106
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('ZET') and not target_is_curated:
        name = name.replace('ZETA', 'ZET ').replace('-', ' ').replace('_', ' ').replace('ZET02RET', 'ZET2 RET').replace('ZET1RET', 'ZET1 RET').replace('ZET2RET', 'ZET2 RET').replace('ZETTUC', 'ZET TUC')
        curation = 107
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('ETA') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ').replace('ETAAPS', 'ETA APS').replace('ETAORI', 'ETA ORI')
        curation = 108
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('THE') and not target_is_curated:
        name = name.replace('-', ' ').replace('THETA', 'THE ')
        curation = 109
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('TET') and not target_is_curated: # THETA variants: TETA
        name = name.replace('-', ' ').replace('TET', 'THE ')
        curation = 110
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('IOT') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ').replace('IOTHOR', 'IOT HOR').replace('IOTCAR', 'IOT CAR')
        curation = 110
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('KAP') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ').replace('KAPPA', 'KAP').replace('KAPFOR', 'KAP FOR').replace('KAP01', 'KAP1')
        curation = 111
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('LAM') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ').replace('LAMBDA', 'LAM')
        curation = 112
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('MU') and not target_is_curated:
        name = name.replace('MUARA', 'mu. ARA').replace('MU1SCO', 'mu.01 Sco').replace('MULUP', 'mu. LUP').replace('MUCEN', 'mu. CEN').replace('MU02CRU', 'mu.02 CRU')
        curation = 113
        target_is_curated = True
        ncur = ncur +1
        if name == 'MUS':
            name = "ETA MUS"
            curation = 113
            target_is_curated = True
    if name.startswith('NU') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ').replace('NU', 'nu.')
        curation = 114
        target_is_curated = True
        ncur = ncur +1
    if name.startswith(('XI','KSI')) and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ').replace('XI', 'KSI')
        curation = 115
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('OMI') and not target_is_curated:
        name = name.replace('-', ' ').replace('0', ' ').replace('OMILUP', 'OMI LUP')
        curation = 116
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('PI') and not target_is_curated:
        name = "pi. HYA"
        curation = 117
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('RHO') and not target_is_curated:
        name = name.replace('-', ' ').replace('RHOPAV', 'RHO PAV')
        curation = 118
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('SIG') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ')
        curation = 119
        target_is_curated = True
    if name.startswith('TAU') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ').replace('CETI', 'CET').replace('TAUCET', 'TAU CET').replace('TAUBOO', 'TAU BOO')
        curation = 120
        target_is_curated = True
        ncur = ncur +1
        if name == 'TAU':
            name = "lam TAU"
            curation = 120
            target_is_curated = True
    if name.startswith('UPS') and not target_is_curated:
        name = "UPS SCO"
        curation = 121
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('PHI') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ').replace('PHI2PAV', 'phi02 PAV')
        curation = 122
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('CHI') and not target_is_curated:
        name = name.replace('CHI-VIR_K2III', 'CHI VIR').replace('-', ' ').replace('_', ' ')
        curation = 123
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('PSI') and not target_is_curated:
        name = name.replace('-', ' ')
        curation = 124
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('OME') and not target_is_curated:
        name = name.replace('-', ' ')
        curation = 125
        target_is_curated = True

    # stars
    if original_name.startswith('NGC') and not target_is_curated:
        if original_name.startswith('NGC '):
            name = original_name
            curation = 2002
            target_is_curated = True
            ncur = ncur +1
            note = 'NGC'
        elif original_name[3].isdigit():
            name = "NGC " + original_name[3:7] + " "
            name += ''.join(char for char in original_name[7:] if char.isdigit())
            curation = 2002
            target_is_curated = True
            ncur = ncur +1
            note = 'NGC'
    if original_name.startswith('IC') and not target_is_curated:
        if original_name[2].isdigit():
            name = "IC " + original_name[2:6] + " "
            name += ''.join(char for char in original_name[6:] if char.isdigit())
            curation = 2003
            target_is_curated = True
            ncur = ncur +1
            note = 'IC'
    if name.startswith('TOI') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ')
        curation = 2004
        target_is_curated = True
        note = 'TESS'
        ncur = ncur +1
    if name.startswith('EPIC') and not target_is_curated:
        name = name.replace('-', ' ')
        curation = 5
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('WASP') and not target_is_curated:
        name = name.replace('_ORB', '').replace('_TRA', '')
        curation = 6
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('TYC') and not target_is_curated:
        name = name.replace('_', ' ').replace('TYC-', 'TYC ')
        curation = 7
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('BD') and not target_is_curated:
        name = name.replace('_', ' ').replace('--', '-')
        curation = 8
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('NLTT') and not target_is_curated:
        name = name.replace('-', ' ').replace('NLTT44619_F2V', 'NLTT44619')
        curation = 9
        target_is_curated = True
        ncur = ncur +1
    if re.match("^HD[0-9]", name) and "_" in name and not target_is_curated:
        name = re.sub('_.*', '', name)
        curation = 10
        target_is_curated = True
        ncur = ncur +1
    if re.match("^HD[0-9]", name) and "-" in name and not target_is_curated:
        name = re.sub('-.*', '', name)
        curation = 11
        target_is_curated = True
        ncur = ncur +1
    if re.match("^HR[0-9]", name) and "_" in name and not target_is_curated:
        name = re.sub('_.*', '', name)
        curation = 12
        target_is_curated = True
        ncur = ncur +1
    if re.match("^HR[0-9]", name) and "-" in name and not target_is_curated:
        name = re.sub('-.*', '', name)
        curation = 13
        target_is_curated = True
        ncur = ncur +1
    if re.match("^HIP[0-9]", name) and "_" in name and not target_is_curated:
        name = re.sub('_.*', '', name)
        curation = 14
        target_is_curated = True
        ncur = ncur +1
    if re.match("^HIP[0-9]", name) and "-" in name and not target_is_curated:
        name = re.sub('-.*', '', name)
        curation = 15
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('V ') and not target_is_curated:
        name = name.replace('V ', 'V* ')
        curation = 16
        target_is_curated = True
        ncur = ncur +1
    if re.match("^V[0-9]", name) and any(name.endswith(c) for c in constellations) and not target_is_curated:
        name = "V* "+name
        name = name.replace('-', ' ').replace('_', ' ')
        pref=name[:-3]
        suff=name[-3:]
        name=pref+" "+suff
        curation = 17
        target_is_curated = True
        ncur = ncur +1
    if re.match(r"(V\d+)([A-Z]{3})", name) and any(name.endswith(c) for c in constellations) and not target_is_curated:
        match = re.match(r"(V\d+)([A-Z]{3})", name)
        name = " ".join(match.groups())
        curation = 18
        target_is_curated = True
        ncur = ncur +1
    if name.startswith(('LRC0', 'LRA0', 'IRA0', 'SRC0', 'SRA0')) and not target_is_curated:
        name = name.replace('_', ' ')
        curation = 2019
        target_is_curated = True
        ncur = ncur +1
        note = "EXODAT"
    if name.startswith('TIC') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ')
        curation = 20
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('HIP') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ')
        suff = name[-1:] 
        if suff.isalpha():
            name = name[:-1]
        curation = 21
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('HD') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ')
        suff = name[-1:] 
        if suff.isalpha() and suff in ['A', 'B', 'C', 'D', 'E', 'F']:
            name = name[:-1]+" "+suff
        elif suff.isalpha() and suff not in ['A', 'B', 'C', 'D', 'E', 'F']:
            name = name[:-1]
        curation = 22
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('HR') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ')
        suff = name[-1:] 
        if suff.isalpha() and suff in ['A', 'B', 'C', 'D', 'E', 'F']:
            name = name[:-1]+" "+suff
        elif suff.isalpha() and suff not in ['A', 'B', 'C', 'D', 'E', 'F']:
            name = name[:-1]
        curation = 23
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('LHS') and not target_is_curated:
        name = name.replace('LHS', 'LHS ')
        if '_' in name:
            name = name[:name.index('_')]
            curation = 24
            target_is_curated = True
            ncur = ncur +1
    if name.startswith('CD') and not target_is_curated:
        name = name.replace('--', '-').replace('CD-53-251', 'CD-53 251').replace('CD-69 1055', 'CD-69 1055').replace('CD-71-1234', 'CD-71 1234').replace('CD-84-0080', 'CD-84 0080').replace('CD-TAU', 'CD TAU').replace('CD_CIR', 'CD CIR')
        curation = 25
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('HV') and not target_is_curated:
        name = name.replace('HV-LUP', 'HV LUP').replace('HV_LUP', 'HV LUP').replace('HVERI', 'HV ERI')
        curation = 26
        target_is_curated = True
        ncur = ncur +1
    if re.fullmatch(r"E_\d{9}", name) and not target_is_curated:
        digits = re.search(r"\d{9}", name).group()
        name = "EPIC " + digits
        curation = 27
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('BLG') and not target_is_curated:
        parts = name.split('_')
        name = "OGLE " + " ".join(parts)
        curation = 2028
        target_is_curated = True
        ncur = ncur +1
        note = 'OGLE'
    if re.fullmatch(r"LMC-CEP-\d{4}", name) and not target_is_curated:
        name = "OGLE " + name.replace('-', ' ')
        curation = 2029
        target_is_curated = True
        ncur = ncur +1
        note = 'OGLE'
    if re.fullmatch(r"LMC\d{3}\.\d-\d+", name) and not target_is_curated:
        parts = name.split('-')
        name = "OGLE " + parts[0] + parts[1]
        curation = 2030
        target_is_curated = True
        ncur = ncur +1
        note = 'OGLE'
    if re.fullmatch(r"LMC\d{3}\.\d{2}[-_.]\d+", name) and not target_is_curated:
        parts = re.split('[-_.]', name)
        name = "OGLE " + parts[0] + " " + parts[1]
        curation = 2031
        target_is_curated = True
        ncur = ncur +1
        note = 'OGLE'
    if name.startswith('OGLE') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ')
        curation = 2032
        target_is_curated = True
        ncur = ncur +1
        note = 'OGLE'
    if re.fullmatch(r"SMC-CEP-\d{4}", name) and not target_is_curated:
        name = "OGLE " + name.replace('-', ' ')
        curation = 2033
        target_is_curated = True
        ncur = ncur +1
        note = 'OGLE'
    if name.startswith(('GL', 'GJ')) and not target_is_curated:
        name = name.replace('-', '').replace('_', '').replace(' ', '')
        curation = 34
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('CPD') and not target_is_curated:
        name = name.replace('+', '-').replace('_', ' ').replace('-', 'X',1).replace('-', ' ',1).replace('X', '-',1)
        curation = 35
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('COROT') and not target_is_curated:
        name = original_name
        curation = 36
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('GSC') and not target_is_curated:
        name = name.replace('-', ' ', 1).replace('_', ' ')
        curation = 37
        target_is_curated = True
        ncur = ncur +1
    if re.fullmatch(r"HATS\d{3}\-\d{3}", name) and not target_is_curated:
        parts = re.split('[-_.]', name)
        name = name[:4]+" "+name[4:]
        note = "HATS"
        curation = 2038
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('LTT') and not target_is_curated:
        name = name.replace('-', ' ')
        curation = 39
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('LHS') and not target_is_curated:
        name = name.replace('-', ' ')
        curation = 40
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('LP') and not target_is_curated:
        name = name.replace('_', ' ')
        curation = 41
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('LMC') and not target_is_curated:
        name = "OGLE " + name.replace('-', ' ')
        curation = 42
        target_is_curated = True
        ncur = ncur +1
        note = 'OGLE'
    if name.startswith('SMC') and not target_is_curated:
        name = "OGLE " + name.replace('-', ' ')
        curation = 43
        target_is_curated = True
        ncur = ncur +1
        note = 'OGLE'
    if name.startswith('UCAC') and not target_is_curated:
        name = name.replace('_', ' ')
        curation = 44
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('V*') and not target_is_curated:
        name = name.replace('_', ' ')
        curation = 45
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('WOLF') and not target_is_curated:
        name = name.replace('-', ' ')
        curation = 46
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('2MASS') and not target_is_curated:
        name = name.replace('_', ' ')
        curation = 47
        target_is_curated = True
        ncur = ncur +1
    if len(name)>=5 and re.match("^[A-Z][A-Z]", name) and any(name.endswith(c) for c in constellations) and not target_is_curated:
        pref=name[:2]
        suff=name[-3:]
        name="V* "+pref+" "+suff
        #print(original_name, " | ", name)
        curation = 48
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('CCDM') and not target_is_curated:
        name = name.replace('-', ' ',1)
        curation = 49
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('M67') and not target_is_curated:
        name = name.replace('-', ' ',1)
        curation = 50
        target_is_curated = True
        ncur = ncur +1
        note = 'M67'
    if name.startswith('MELOTTE71') and not target_is_curated:
        name = name.replace('NO.', ' ',1)
        curation = 51
        target_is_curated = True
        ncur = ncur +1
        note = 'MELOTTE71'
    if len(name)>=4 and re.match("^[A-Z]", name) and any(name.endswith(c) for c in constellations) and not target_is_curated:
        pref=name[:1]
        suff=name[-3:]
        name="V* "+pref+" "+suff
        #print(original_name, " | ", name)
        curation = 52
        target_is_curated = True
        ncur = ncur +1
    if re.match("^[0-9][0-9][0-9][0-9]_[A-Z][A-Z][-_]", name) and any(name.endswith(c) for c in constellations) and not target_is_curated:
        pref=name[5:7]
        suff=name[-3:]
        name='V* '+pref+" "+suff
        #print('aa ',original_name, " | ", name)
        curation = 53
        target_is_curated = True
        ncur = ncur +1
    if re.match("^[0-9][0-9][0-9][0-9]_V[0-9][0-9][0-9][-_]", name) and any(name.endswith(c) for c in constellations) and not target_is_curated:
        pref=name[5:9]
        suff=name[-3:]
        name='V* '+pref+" "+suff
        #print('bb ',original_name, " | ", name)
        curation = 53
        target_is_curated = True
        ncur = ncur +1
    if re.match("^[0-9][0-9][0-9][0-9]_V[0-9][0-9][0-9][0-9][-_]", name) and any(name.endswith(c) for c in constellations) and not target_is_curated:
        pref=name[5:10]
        suff=name[-3:]
        name='V* '+pref+" "+suff
        #print('cc ',original_name, " | ", name)
        curation = 53
        target_is_curated = True
        ncur = ncur +1
    if len(name)>=6 and re.match("^[0-9][0-9][0-9]", name) and any(name.endswith(c) for c in constellations) and not target_is_curated:
        pref=name[:3]
        suff=name[-3:]
        name=pref+" "+suff
        #print('x1',original_name, " | ", name)
        curation = 54
        target_is_curated = True
        ncur = ncur +1
    if len(name)>=5 and re.match("^[0-9][0-9]", name) and any(name.endswith(c) for c in constellations) and not target_is_curated:
        pref=name[:2]
        suff=name[-3:]
        name=pref+" "+suff
        #print('x2',original_name, " | ", name)
        curation = 55
        target_is_curated = True
        ncur = ncur +1
    if len(name)>=4 and re.match("^[0-9]", name) and any(name.endswith(c) for c in constellations) and not target_is_curated:
        pref=name[:1]
        suff=name[-3:]
        name=pref+" "+suff
        #print('x3',original_name, " | ", name)
        curation = 56
        target_is_curated = True
        ncur = ncur +1
    if len(name)>=4 and re.match("^[A-Z]", name) and any(name.endswith(c) for c in constellations) and not target_is_curated:
        pref=name[:1]
        suff=name[-3:]
        name=pref+" "+suff
        #print('yy',original_name, " | ", name)
        curation = 57
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('KELT') and not target_is_curated:
        name = name.replace('_', '-')
        curation = 58
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('HAT') and not target_is_curated:
        name = name.replace('-', ' ',1).replace('_', ' ')
        curation = 59
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('TWA') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ')
        curation = 60
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('SAO') and not target_is_curated:
        name = name.replace('-', ' ').replace('_', ' ')
        curation = 61
        target_is_curated = True
        ncur = ncur +1
    if name.startswith(('QSO','Q0515','HE-0515','HE0515')) and not target_is_curated:
        name = 'QSO B0515-4414'
        curation = 62
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('K2') and not target_is_curated:
        name = original_name
        curation = 63
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('IRAS') and not target_is_curated:
        name = name.replace('_', ' ')
        curation = 64
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('G') and not target_is_curated:
        nome = name.replace(' ', '')
        if re.match("^G[0-9][0-9][0-9]-[0-9]", nome):
            name=nome
            curation = 65
            target_is_curated = True
            ncur = ncur +1
    if name.startswith('L') and not target_is_curated:
        nome = name.replace(' ', '')
        if re.match("^L[0-9][0-9][0-9]-[0-9]", nome):
            name=nome
            curation = 66
            target_is_curated = True
            ncur = ncur +1



    if name.startswith('J') and not target_is_curated:
        name = original_name
        note = "2MASS"
        curation = 1001
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('SERAM') and not target_is_curated:
        name = original_name
        note = "2MASS"
        curation = 1002
        target_is_curated = True
        ncur = ncur +1
    if name.startswith('SW') and name[2].isdigit():
        name = original_name
        note = "2MASS"
        curation = 1003
        target_is_curated = True
        ncur = ncur +1



    # Replace characters as needed
    #name = name.replace('_', ' ', 1).replace('--', ' ')
    # Collapse multiple spaces into one
    name = ' '.join(name.split())
    if curation == 0 and kuration == 0:
        name = "nc"
        note = ""
    if ncur>1:
        print(ncur, " | ", original_name, " | ", name)

    return name, curation,kuration,note


# Create lists to store the results
target_curated = []
curation = []
curation_basic = []
note = []


# Iterate over the rows of the DataFrame
for i in range(len(df)):
    # Apply the function to the first column
    name=df.iloc[i, 0]
    target_curated_, curation_, kuration_, notes_ = curate_name(name)
    # Append the results to the lists
    target_curated.append(target_curated_)
    curation.append(curation_)
    curation_basic.append(kuration_)
    note.append(notes_)

# Insert the results as the second and third columns
df.insert(1, 'target_curated', target_curated)
df.insert(2, 'curation', curation)
df.insert(3, 'curation_basic', curation_basic)
df.insert(4, 'note', note)

cols = df.columns.tolist()
idx_e_ra = cols.index('e_ra')
idx_dec = cols.index('dec')
cols[idx_e_ra], cols[idx_dec] = cols[idx_dec], cols[idx_e_ra]
df = df[cols]

# Sort the DataFrame by the second column
df.sort_values(by=df.columns[1], inplace=True)

#df_filtered = df[df['curation'] == 0]
#df_filtered.to_csv('targets1.csv', index=False)


df_filtered = df[(df['curation'] > 0) & (df['curation'] < 1000) & (df['note'] == "")]
df_filtered_unique = df_filtered['target_curated'].drop_duplicates()
df_filtered_unique.to_csv('curation_output/targets_for_simbad.csv', index=False, header=False)
df_filtered_unique = df_filtered.drop_duplicates(subset='target_curated')
df_filtered_unique.to_csv('curation_output/targets_for_simbad.all.csv', index=False, header=True)


df_filtered = df[(df['curation'] > 1000) & (df['curation'] < 10000)]
df_filtered_unique = df_filtered['target_curated'].drop_duplicates()
df_filtered_unique.to_csv('curation_output/targets_for_vizier.csv', index=False, header=False)
df_filtered_unique = df_filtered.drop_duplicates(subset='target_curated')
df_filtered_unique.to_csv('curation_output/targets_for_vizier.all.csv', index=False, header=True)

df_filtered = df[df['curation'] == 0]
df_filtered_unique = df_filtered.drop_duplicates(subset='target')
df_filtered_unique.to_csv('curation_output/targets_not_simbad.all.csv', index=False, header=True)

df_filtered = df[(df['curation'] > 100000)]
df_filtered_unique = df_filtered.drop_duplicates(subset='target_curated')
df_filtered_unique.to_csv('curation_output/targets_for_horizons.all.csv', index=False, header=True)


# Write the result to a new CSV file
#df.to_csv('targets1.csv', index=False)

nt=df.shape[0]
nc=df['curation'].value_counts().get(0, 0)
nm=nt-nc
print("TARGETS    ",nt)
print("CURATED    ",nm)
print("NOT CURATED",nc)
