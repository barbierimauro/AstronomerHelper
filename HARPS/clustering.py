#!/scratch/mbarbier/miniconda3/bin/python3
import pandas as pd
from sklearn.cluster import DBSCAN
import numpy as np
from scipy.stats import kurtosis, skew
import sys
from numpy.linalg import eig
from shapely.geometry import MultiPoint, Point
from shapely.ops import cascaded_union
from astropy.time import Time
from collections import Counter
import warnings

# Ignore all warnings
warnings.filterwarnings('ignore')

def top_n_labels(series,n=5):
    # Get a list of the counts of each label
    counts = Counter(series).most_common(n)
    
    # If there are less than n labels, extend the counts list with empty strings
    while len(counts) < n:
        counts.append(("", 0))
        
    # Handle "", "NONE", and "nan" labels
    cleaned_counts = [(str(label) if label not in ["", "NONE", "nan"] else "", _) for label, _ in counts]
    
    # Prepare the string for return
    result = ';'.join([label for label, _ in cleaned_counts])
    
    return result

def count_ordered(group):

    order_dp_tech = ["ECHELLE","ECHELLE,CIRPOL","ECHELLE,LINPOL","ECHELLE,ABSORPTION-CELL"]
    order_wl_res_power = [115000, 80000]
    order_dp_cat = ["SCIENCE","CALIB"]

    # Get the counts for 'dp_tech'
    dp_tech_counts = [group['dp_tech'].value_counts().get(key, 0) for key in order_dp_tech]

    # Count the non-null values in 'archive_id_spec'
    archive_id_spec_count = [group['archive_id_spectra'].notnull().sum()]

    # Get the counts for 'wl_res_power'
    wl_res_power_counts = [group['wl_res_power'].value_counts().get(key, 0) for key in order_wl_res_power]

    # Get the counts for 'dp_cat'
    dp_cat_counts = [group['dp_cat'].value_counts().get(key, 0) for key in order_dp_cat]

    # Combine all counts into a list
    all_counts = dp_tech_counts + archive_id_spec_count + wl_res_power_counts + dp_cat_counts

    # Define column names
    column_names = ['EC','CP','LP','I2','P3','HAM','EGG', 'SCI', 'CAL']

    # Return the counts as a series with the specified column names
    return pd.Series(all_counts, index=column_names)

# function to compute the bounding box and its spatial density
def compute_bbox(group):
    points = [Point(row['ra'], row['dec']) for _, row in group.iterrows()]
    multipoint = MultiPoint(points)
    bbox = multipoint.convex_hull.envelope.exterior.coords[:-1]  # compute the bounding box corners
    bbox_string = ";".join([f"{x:.6f};{y:.6f}" for (x, y) in bbox])
    # Calculate width, height and diagonal of the bounding box
    bbox_width = (bbox[2][0] - bbox[0][0])*3600
    bbox_height = (bbox[2][1] - bbox[0][1])*3600
    bbox_diagonal = np.sqrt(bbox_width**2 + bbox_height**2)
    bbox_area = bbox_width * bbox_height
    if bbox_area > 0:
        spatial_density = len(points) / bbox_area
    else:
        spatial_density = np.nan
    
    bbox_output = {
        'x1': bbox[0][0],
        'y1': bbox[0][1],
        'x2': bbox[1][0],
        'y2': bbox[1][1],
        'x3': bbox[2][0],
        'y3': bbox[2][1],
        'x4': bbox[3][0],
        'y4': bbox[3][1],
        'bbox_width': bbox_width, 
        'bbox_height': bbox_height, 
        'bbox_diagonal': bbox_diagonal, 
        'spatial_density': spatial_density
    }
       
    return pd.Series(bbox_output)


# function to compute the ellipse parameters
def compute_ellipse_params(group):
    if len(group) > 3:
        # compute the covariance matrix
        cov = np.cov(group['ra'], group['dec'])
        # compute the eigenvalues and eigenvectors
        eigenvalues, eigenvectors = eig(cov)
        # compute the lengths of the semi-major and semi-minor axes
        a = np.sqrt(5.991 * max(eigenvalues))
        b = np.sqrt(5.991 * min(eigenvalues))
        # compute the orientation
        angle = np.arctan2(*eigenvectors[:, np.argmax(eigenvalues)]) * 180 / np.pi
    else:
        a, b, angle = np.nan, np.nan, np.nan
    return pd.Series({
        'semi_major_axis': a,
        'semi_minor_axis': b,
        'position_angle': angle
    })

def clustering(dataframe,filename,epsilon):

    filestat=filename+".stats.csv"
    filecat =filename+".cat.csv"
    points = dataframe[['ra', 'dec']]

    # perform DBSCAN clustering
    eps_dbscan = epsilon
    dbscan = DBSCAN(eps=eps_dbscan, min_samples=2)
    dbscan_labels = dbscan.fit_predict(points)

    # save DBSCAN results to a new DataFrame
    dbscan_results = dataframe.assign(DBSCAN_Label=dbscan_labels)
    
    
    # calculate statistics for each cluster
    dbscan_stats = dbscan_results.groupby('DBSCAN_Label').agg({
        'ra': ['mean', 'median', 'std', skew, kurtosis, 'count'],
        'dec': ['mean', 'median', 'std', skew, kurtosis]
    })


    mode_object = dbscan_results.groupby('DBSCAN_Label')['object'].apply(top_n_labels,n=11).reset_index()
    mode_object.columns = ['DBSCAN_Label', 'object']
    
    mode_spt_obj = dbscan_results.groupby('DBSCAN_Label')['spt_obj'].apply(top_n_labels,n=11).reset_index()
    mode_spt_obj.columns = ['DBSCAN_Label', 'spt']
    
    mode_title = dbscan_results.groupby('DBSCAN_Label')['prog_title'].apply(top_n_labels,n=1).reset_index()
    mode_title.columns = ['DBSCAN_Label', 'prog_title']
    
    mode_progid = dbscan_results.groupby('DBSCAN_Label')['prog_id'].apply(top_n_labels,n=1).reset_index()
    mode_progid.columns = ['DBSCAN_Label', 'prog_id']
    
    dbscan_results['picoi1'] = dbscan_results['pi_coi'].str.split('/').str[0]
    mode_picoi = dbscan_results.groupby('DBSCAN_Label')['picoi1'].apply(top_n_labels,n=1).reset_index()
    mode_picoi.columns = ['DBSCAN_Label', 'picoi1']
    
    
    # compute the ellipse parameters for each cluster
    ellipse_params = dbscan_results.groupby('DBSCAN_Label').apply(compute_ellipse_params).reset_index()
    
    # compute the bounding box and spatial density for each cluster
    bbox = dbscan_results.groupby('DBSCAN_Label').apply(compute_bbox).reset_index()
    
    # compute the closest point to the median for each cluster
    dbscan_results['dist_to_median'] = np.sqrt((dbscan_results['ra'] - dbscan_results['ra'].median())**2 + (dbscan_results['dec'] - dbscan_results['dec'].median())**2)
    median_mjd = dbscan_results.loc[dbscan_results.groupby('DBSCAN_Label')['dist_to_median'].idxmin()][['DBSCAN_Label', 'mjd_obs']]
    median_mjd.columns = ['DBSCAN_Label', 'median_mjd']
    
    # compute pseudomag for each cluster
    dbscan_results['pf'] = np.where((dbscan_results['snr'] > 0) & (dbscan_results['snr'] < 1000) & (dbscan_results['exptime'] > 0) & (dbscan_results['exptime'] < 40000), dbscan_results['snr'] / dbscan_results['exptime'], np.nan)
    pseudomag = dbscan_results.groupby('DBSCAN_Label')['pf'].median().reset_index()
    pseudomag['pseudomag'] = -3.3 * np.log10(pseudomag['pf']) + 5.5
    pseudomag = pseudomag[['DBSCAN_Label', 'pseudomag']]
    
    # time
    # min
    time_min = dbscan_results.groupby('DBSCAN_Label')['mjd_obs'].min().reset_index()
    time_min.columns = ['DBSCAN_Label', 'min_mjd']
    # max
    time_max = dbscan_results.groupby('DBSCAN_Label')['mjd_obs'].max().reset_index()
    time_max.columns = ['DBSCAN_Label', 'max_mjd']
    # range
    time_span = dbscan_results.groupby('DBSCAN_Label')['mjd_obs'].apply(lambda x: x.max() - x.min()).reset_index()
    time_span.columns = ['DBSCAN_Label', 'time_span']
    
    #counts
    obs_counts = dbscan_results.groupby('DBSCAN_Label').apply(count_ordered).reset_index()
    
    
    # flatten the multi-level column index
    dbscan_stats.columns = ['_'.join(col) for col in dbscan_stats.columns]
    
    # merge the computed parameters with the statistics DataFrame
    dbscan_stats = dbscan_stats.reset_index().join(ellipse_params.set_index('DBSCAN_Label'), on='DBSCAN_Label')
    dbscan_stats = dbscan_stats.join(bbox.set_index('DBSCAN_Label'), on='DBSCAN_Label')
    dbscan_stats = dbscan_stats.merge(median_mjd,   left_on='DBSCAN_Label', right_on='DBSCAN_Label')
    dbscan_stats = dbscan_stats.merge(pseudomag,    left_on='DBSCAN_Label', right_on='DBSCAN_Label')
    dbscan_stats = dbscan_stats.merge(mode_object,  left_on='DBSCAN_Label', right_on='DBSCAN_Label', how='left')
    dbscan_stats = dbscan_stats.merge(mode_spt_obj, left_on='DBSCAN_Label', right_on='DBSCAN_Label', how='left')
    dbscan_stats = dbscan_stats.merge(time_min,     left_on='DBSCAN_Label', right_on='DBSCAN_Label', how='left')
    dbscan_stats = dbscan_stats.merge(time_max,     left_on='DBSCAN_Label', right_on='DBSCAN_Label', how='left')
    dbscan_stats = dbscan_stats.merge(time_span,    left_on='DBSCAN_Label', right_on='DBSCAN_Label', how='left') 
    dbscan_stats = dbscan_stats.merge(mode_title,   left_on='DBSCAN_Label', right_on='DBSCAN_Label', how='left') 
    dbscan_stats = dbscan_stats.merge(mode_progid,  left_on='DBSCAN_Label', right_on='DBSCAN_Label', how='left') 
    dbscan_stats = dbscan_stats.merge(mode_picoi,   left_on='DBSCAN_Label', right_on='DBSCAN_Label', how='left') 
    dbscan_stats = dbscan_stats.merge(obs_counts,   left_on='DBSCAN_Label', right_on='DBSCAN_Label', how='left')
    
    
    # multiply the errors by 3600
    dbscan_stats['ra_std']  *= 3600
    dbscan_stats['dec_std'] *= 3600
    # multiply the semi-major and semi-minor axes by 3600
    dbscan_stats['semi_major_axis'] *= 3600
    dbscan_stats['semi_minor_axis'] *= 3600
    
    # rename columns
    dbscan_stats.rename(columns={
        'DBSCAN_Label': 'GroupID',
        'ra_count': 'npoints'
    }, inplace=True)
    
    # Remove the row with GroupID == -1
    
    dbscan_stats = dbscan_stats[dbscan_stats.GroupID != -1]
    
    
    # order the columns
    dbscan_stats = dbscan_stats[[
        'GroupID', 'object', 'spt',
        'npoints', 'median_mjd', 'pseudomag',
        'ra_mean', 'dec_mean', 'ra_median', 'dec_median',
        'ra_std', 'dec_std', 'ra_skew', 'dec_skew',
        'ra_kurtosis', 'dec_kurtosis',
        'semi_major_axis', 'semi_minor_axis', 'position_angle',
        'bbox_width', 'bbox_height', 'bbox_diagonal',
        'spatial_density',
        'time_span', 'min_mjd', 'max_mjd',
        'x1','y1','x2','y2','x3','y3','x4','y4',
        'EC','CP','LP','I2','P3','HAM','EGG', 'SCI', 'CAL',
        'prog_id', 'picoi1', 'prog_title'
    ]]

    groups_num=len(dbscan_stats)
    grouped_total = dbscan_stats['npoints'].sum()
    total_objects=len(points)
    num_removed=total_objects-grouped_total
    print()
    print('subset       :',filename)
    print('N groups     :',groups_num)
    print('total grouped:',grouped_total)
    print('not grouped  :',num_removed)
    print()

    
    # save DBSCAN statistics to a new CSV file
    dbscan_stats.to_csv(filestat, index=False)

    # save only the unclustered results
    unclustered_results = dbscan_results.loc[dbscan_results['DBSCAN_Label'] == -1]

    # keep only certain columns
    columns_to_keep = ['dp_id_raw','ob_name','object','target','ra','dec','pseudo_hpmag']
    unclustered_results = unclustered_results[columns_to_keep]

    # Then write to csv
    unclustered_results.to_csv(filecat, index=False)



    return

########################################################################
########################################################################
########################################################################
########################################################################





# load the data
df = pd.read_csv('FINAL_all_observations_harps.fits.csv')
df['object'] = df['object'].fillna('')
df['prog_title'] = df['prog_title'].fillna('')

#print(len(df))
#sys.exit()
#df['object_original'] = df['object']

df.loc[df['object'].str.startswith('STAR'), 'object'] = df['target']
df.loc[df['ob_name'].str.startswith('Wasp17'), 'object'] = df['ob_name']
df.loc[df['ob_name'].str.startswith('WCEN'), 'object'] = 'Ome Cen'
df.loc[(df['ob_name'].str.strip() == '') & (df['target'].str.startswith('HD') | df['target'].str.startswith('HIP')), 'object'] = df['target']
df.loc[df['object'].isin(['NO NAME', 'NO-NAME']), 'object'] = df['ob_name']

df.loc[df['object'].str.lower().str.contains('sun'), 'object']       = 'SUN'
df.loc[df['object'].str.lower().str.contains('solar'), 'object']     = 'SUN'
df.loc[df['object'].str.lower().str.contains('sky'), 'object']       = 'SUN'
df.loc[df['object'].str.lower().str.contains('moon'), 'object']      = 'MOON'
df.loc[df['object'].str.lower().str.contains('mooon'), 'object']     = 'MOON'
df.loc[df['object'].str.lower().str.contains('tycho'), 'object']     = 'MOON'
df.loc[df['object'].str.lower().str.contains('giove'), 'object']     = 'JUPITER'
df.loc[df['object'].str.lower().str.contains('jupiter'), 'object']   = 'JUPITER'
df.loc[df['object'].str.lower().str.contains('ganymede'), 'object']  = 'GANYMEDE'
df.loc[df['object'].str.lower().str.contains('ganimede'), 'object']  = 'GANYMEDE'
df.loc[df['object'].str.lower().str.contains('6-hebe'), 'object']    = 'HEBE'
df.loc[df['object'].str.lower().str.contains('comet'), 'object']     = '73P-C/Schwassmann-Wachmann 3'
df.loc[df['object'].str.lower().str.contains('ceres'), 'object']     = 'CERES'
df.loc[df['object'].str.lower().str.contains('vesta'), 'object']     = 'VESTA'
df.loc[df['object'].str.lower().str.contains('astraea'), 'object']   = 'ASTRAEA'
df.loc[df['object'].str.lower().str.contains('europa'), 'object']    = 'EUROPA'
df.loc[df['object'].str.lower().str.contains('herculina'), 'object'] = 'HERCULINA'
df.loc[df['object'].str.lower().str.contains('irene'), 'object']     = 'IRENE'
df.loc[df['object'].str.lower().str.contains('iris'), 'object']      = 'IRIS'
df.loc[df['object'].str.lower().str.contains('juno'), 'object']      = 'JUNO'
df.loc[df['object'].str.lower().str.contains('massalia'), 'object']  = 'MASSALIA'
df.loc[df['object'].str.lower().str.contains('melpomene'), 'object'] = 'MELPOMENE'
df.loc[df['object'].str.lower().str.contains('pallas'), 'object']    = 'PALLAS'
df.loc[df['object'].str.lower().str.contains('venus'), 'object']     = 'VENUS'
df.loc[df['object'].str.lower().str.contains('titan'), 'object']     = 'TITAN'
df.loc[df['object'].str.lower() == 'io', 'object']                   = 'IO'

df.loc[df['object'].str.upper().str.contains('HD10700|TAU-CET|HR509|TAU CET|TAUCETI|TAU-CETI|TAUCET|HD-10700|HD 10700'), 'object']     = 'TAU CET'
df.loc[df['object'].str.upper().str.contains('HD039060|HD39060|BETA-PIC|BETA PIC|HIP 27321|BETA_PIC|BET-PIC'), 'object']     = 'BET PIC'
df.loc[df['object'].str.upper().str.contains('HD20794|E ERI'), 'object']     = 'HD 20794'
df.loc[df['object'].str.upper().str.contains('HD190248|HD-190248|DELTA-PAV|DEL-PAV|HR-7665'), 'object']     = 'DEL PAV'
df.loc[df['object'].str.upper().str.contains('18_SCO|18 SCO|HD146233|18SCO|HIP-79672|HIP79672|HD 146233'), 'object']     = '18 SCO'
df.loc[df['object'].str.upper().str.contains('PROCYON|ALF_CMI'), 'object']     = 'ALF CMI'
df.loc[df['object'].str.upper().str.contains('ALP-CIR|HD128898|HD-128898-V|HD-128898'), 'object']     = 'ALF CIR'
df.loc[df['object'].str.upper().str.contains('HD121370;ETA-BOO'), 'object']     = 'ETA BOO'
df.loc[df['object'].str.upper().str.contains('EPS-INDI-A;HD209100;EPS IND;HR8387;EPSIND'), 'object']     = 'EPS IND'
df.loc[df['object'].str.upper().str.contains('HD1581;ZET TUC;HR77;ZETTUC'), 'object']     = 'ZET TUC'
df.loc[df['object'].str.upper().str.contains('HD17051|IOT-HOR|HR810|IOTHOR'), 'object']     = 'IOT HOR'
df.loc[df['object'].str.upper().str.contains('BETA-HYI|HR98|HD2151|BET-HYI|HD-2151|BETHYI'), 'object']     = 'BET HYI'
df.loc[df['object'].str.upper().str.contains('HD160691|HD-160691|MUARA'), 'object']     = 'MU. ARA'
df.loc[df['object'].str.upper().str.contains('HD188512|BETA-AQL'), 'object']     = 'BET AQL'
df.loc[df['object'].str.upper().str.contains('HD19994|HD 19994|94 CET'), 'object']     = '94 CET'
df.loc[df['object'].str.upper().str.contains('HD-209625|HD209625'), 'object']     = 'HD 209625'
df.loc[df['object'].str.upper().str.contains('HD216956|ALPHA-PSA|ALF-PSA'), 'object']     = 'ALF PSA'
df.loc[df['object'].str.upper().str.contains('HD85512|HD 85512|HD  85512'), 'object']     = 'HD 85512'
df.loc[df['object'].str.upper().str.contains('HD115617|61VIR|61 VIR|LHS349_G5V'), 'object']     = '61 VIR'
df.loc[df['object'].str.upper().str.contains('HD49933|HD049933'), 'object']     = 'HD 49933'
df.loc[df['object'].str.upper().str.contains('HD23249|HR1136|DELERI'), 'object']     = 'DEL ERI'
df.loc[df['object'].str.upper().str.contains('HD 199288|HD199288'), 'object']     = 'HD 199288'
df.loc[df['object'].str.upper().str.contains('HD48915|SIRIUS'), 'object']     = 'ALF CMA'
df.loc[df['object'].str.upper().str.contains('HD46375|HD46735'), 'object']     = 'HD 46375'
df.loc[df['object'].str.upper().str.contains('HD-101065|HD101065_8_0|PRZYBYLSKIS|HD101065|HD_101065_F8_8.0'), 'object']     = 'HD 101065'
df.loc[df['object'].str.upper().str.contains('BETCEN|HD-122451|BET-CEN'), 'object']     = 'BET CEN'
df.loc[df['object'].str.upper().str.contains('HD65907A|HD65907|HD  65907|HR3454'), 'object']     = 'HD 65907'
df.loc[df['object'].str.upper().str.contains('HD22049|EPSILON ERIDANI|HD-22049|HR1084|EPSILON-ERIDANI|EPSERI'), 'object']     = 'EPS ERI'
df.loc[df['object'].str.upper().str.contains('HD26965|HD26965A|OMI02 ERI|GL166C'), 'object']     = 'HIP 19849'
df.loc[df['object'].str.upper().str.contains('GJ887|GJ-887|HD217987|GJ 887|HD 217987|GL887'), 'object']     = 'HD217987'
df.loc[df['object'].str.upper().str.contains('HD69830|HR3259'), 'object']     = 'HD69830'
df.loc[df['object'].str.upper().str.contains('HD114613|HR4979'), 'object']     = 'HD 114613'
df.loc[df['object'].str.upper().str.contains('HD72673|HD 72673|HD72673_STD'), 'object']     = 'HD72673'
df.loc[df['object'].str.upper().str.contains('ALPHADOR|ALF-DOR|HD29305'), 'object']     = 'ALF DOR'
df.loc[df['object'].str.upper().str.contains('HD207129|HR8323|HD 207129'), 'object']     = 'HD 207129'
df.loc[df['object'].str.upper().str.contains('HD59468|HD 59468'), 'object']     = 'hd 59468'
df.loc[df['object'].str.upper().str.contains('PROXIMA|GL551|PROXIMA-CENTAURI|NAME PROXIMA CENTAURI|GJ551|PROXIMA CENTAURI'), 'object']     = 'PROXIMA CENTAURI'
df.loc[df['object'].str.upper().str.contains('HD1461|HD   1461'), 'object']     = 'HD 1461'
df.loc[df['object'].str.upper().str.contains('GJ588|GJ-588|GL588'), 'object']     = 'HIP 76074 '
df.loc[df['object'].str.upper().str.contains('HD67523|HR-3185'), 'object']     = 'HD 67523'
df.loc[df['object'].str.upper().str.contains('GJ367;GL367'), 'object']     = 'HIP 47780'
df.loc[df['object'].str.upper().str.contains('WASP-189;SW1502-0301'), 'object']     = 'HIP 73608'
df.loc[df['object'].str.upper().str.contains('HD20807;ZET02RET;HR1010;HD 20807;ZET2RET'), 'object']     = 'HIP 15371'
df.loc[df['object'].str.upper().str.contains('GJ9827;C12_9858;BD-025958'), 'object']     = 'HIP 115752'
df.loc[df['object'].str.upper().str.contains('HD 33793;GL191;GJ-191;GJ 191'), 'object']     = 'HD 33793'
df.loc[df['object'].str.upper().str.contains('GL674;GJ-674;GJ674;GJ 674'), 'object']     = 'HIP 85523'
df.loc[df['object'].str.upper().str.contains('HD45184;HD  45184;HD-45184;HIP30503'), 'object']     = 'HD45184'
df.loc[df['object'].str.upper().str.contains('EPIC211311380'), 'object']     = 'HIP 41378'
df.loc[df['object'].str.upper().str.contains('HD 176986;HD176986'), 'object']     = 'HD 176986'
df.loc[df['object'].str.upper().str.contains('HD   4915;HD4915'), 'object']     = 'HD 4915'
df.loc[df['object'].str.upper().str.contains('HD68978A'), 'object']     = 'HD 68978'
df.loc[df['object'].str.upper().str.contains('HD-104237;DX_CHA'), 'object']     = 'HD 104237'
df.loc[df['object'].str.upper().str.contains('GL699;GJ-699;BARNARDS-STAR;NAME BARNARDS STAR;GJ699'), 'object']     = 'BARNARD STAR'
df.loc[df['object'].str.upper().str.contains('KELT11'), 'object']     = 'HD 93396'
df.loc[df['object'].str.upper().str.contains('GAMCEN;CCDM-J12415-4858AB'), 'object']     = 'GAM CEN'
df.loc[df['object'].str.upper().str.contains('ACHERNAR;HD-10144;HR472;ALPERI'), 'object']     = 'ALP ERI'
df.loc[df['object'].str.upper().str.contains('GL54.1;GJ-54.1'), 'object']     = 'HIP 5643'
df.loc[df['object'].str.upper().str.contains('HD-24712;HD24712;HR-1217;HR1217'), 'object']     = 'HD 24712'
df.loc[df['object'].str.upper().str.contains('HD 161098;HD161098'), 'object']     = 'HD161098'
df.loc[df['object'].str.upper().str.contains('GL273;GJ-273'), 'object']     = 'HIP 36208'
df.loc[df['object'].str.upper().str.contains('HD16417;HD-16417'), 'object']     = 'HD 16417'
df.loc[df['object'].str.upper().str.contains('51PEG;HD217014'), 'object']     = '51 PEG'
df.loc[df['object'].str.upper().str.contains('HD 197481;V AU MIC;HD197481;AU-MIC;AUMIC;GL803'), 'object']     = 'HD 197481'
df.loc[df['object'].str.upper().str.contains('GL447;ROSS  128'), 'object']     = 'HIP 57548'
df.loc[df['object'].str.upper().str.contains('HD210918;HIP109821'), 'object']     = 'HD 210918'
df.loc[df['object'].str.upper().str.contains('HD11977;HD-11977'), 'object']     = 'HD 11977'
df.loc[df['object'].str.upper().str.contains('LHS1140'), 'object']     = 'Gaia DR3 2371032916186181760'
df.loc[df['object'].str.upper().str.contains('HD39194;HD 39194'), 'object']     = 'HD 39194'
df.loc[df['object'].str.upper().str.contains('12-BOO;'), 'object']     = '12 BOO'
df.loc[df['object'].str.upper().str.contains('HD45348;CANOPUS;'), 'object']     = 'ALP CAR'
df.loc[df['object'].str.upper().str.contains('HD43834;HR2261;ALPMEN;'), 'object']     = 'ALP MEN'
df.loc[df['object'].str.upper().str.contains('GL581;GJ581;GI581'), 'object']     = 'HIP 74995'
df.loc[df['object'].str.upper().str.contains('SW0146+0242;WASP-76'), 'object']     = 'TYC 32-111-1'
df.loc[df['object'].str.upper().str.contains('HD071155;HD71155'), 'object']     = 'HD 71155'
df.loc[df['object'].str.upper().str.contains('GL667C;GJ667C'), 'object']     = 'Gaia DR3 5975663354131618304'
df.loc[df['object'].str.upper().str.contains('HD38858;HD  38858'), 'object']     = 'HD 38858'
df.loc[df['object'].str.upper().str.contains('GL876'), 'object']     = 'HIP 113020'
df.loc[df['object'].str.upper().str.contains('HD106315;EPIC201437844'), 'object']     = 'HD 106315'
df.loc[df['object'].str.upper().str.contains('HD82943;HD-82943'), 'object']     = 'HD 82943'
df.loc[df['object'].str.upper().str.contains('LQ-HYA;V LQ HYA;HD-82558;HD082558;HD82558'), 'object']     = 'HD 82558'
df.loc[df['object'].str.upper().str.contains('HD 220339;HD220339'), 'object']     = 'HD 220339'
df.loc[df['object'].str.upper().str.contains('EPIC205071984;C02T004'), 'object']     = 'Gaia DR3 4130539180358512768'
df.loc[df['object'].str.upper().str.contains('HD177565;HD 177565'), 'object']     = 'HD 177565'
df.loc[df['object'].str.upper().str.contains('HR-3748;HD-81797;ALF_HYA'), 'object']     = 'ALF HYA'
df.loc[df['object'].str.upper().str.contains('HD16160;HD 16160'), 'object']     = 'HD 16160'
df.loc[df['object'].str.upper().str.contains('HIP-85042;HD157347;HIP85042;HD157347_STD;HD 157347'), 'object']     = 'HD 157347'
df.loc[df['object'].str.upper().str.contains('HD 170493;HD170493'), 'object']     = 'HD 170493'
df.loc[df['object'].str.upper().str.contains('NU_OCT;HR-8254;HD-205478'), 'object']     = 'HD 205478'
df.loc[df['object'].str.upper().str.contains('HD-74874;HD74874'), 'object']     = 'HD 74874'
df.loc[df['object'].str.upper().str.contains('HD-36705;V AB DOR;HD36705'), 'object']     = 'HD 36705'
df.loc[df['object'].str.upper().str.contains('HD181433;HD 181433'), 'object']     = 'HD 181433'
df.loc[df['object'].str.upper().str.contains('GJ3293'), 'object']     = 'Gaia DR3 4893118771316702720'
df.loc[df['object'].str.upper().str.contains('EPIC245950175'), 'object']     = 'Gaia DR2 2413596935442139520'
df.loc[df['object'].str.upper().str.contains('HD22879;HD-22879'), 'object']     = 'HD 22879'
df.loc[df['object'].str.upper().str.contains('GJ163;HIP19394'), 'object']     = 'HIP 19394'
df.loc[df['object'].str.upper().str.contains('HD220507;HIP115577'), 'object']     = 'HIP 115577'
df.loc[df['object'].str.upper().str.contains('HD   8389;HD8389A;HD8389;NLTT4599;'), 'object']     = 'HD 8389'
df.loc[df['object'].str.upper().str.contains('GJ3138;BD-170400'), 'object']     = 'HIP 10037'
df.loc[df['object'].str.upper().str.contains('HD68146;HD068146'), 'object']     = 'HD 68146'
df.loc[df['object'].str.upper().str.contains('HD32564;HD 32564'), 'object']     = 'HD 32564'
df.loc[df['object'].str.upper().str.contains('GL176;GJ176'), 'object']     = 'HD 285968'
df.loc[df['object'].str.upper().str.contains('GL701;HD 165222'), 'object']     = 'HD 165222'
df.loc[df['object'].str.upper().str.contains('GL229'), 'object']     = 'HD 42581'
df.loc[df['object'].str.upper().str.contains('HD10647;HR506'), 'object']     = 'HD 10647'
df.loc[df['object'].str.upper().str.contains('ALPHA_FOR;HR963;ALPFOR'), 'object']     = 'ALF FOR'
df.loc[df['object'].str.upper().str.contains('HD216435;TAU01 GRU'), 'object']     = 'HD 216435'
df.loc[df['object'].str.upper().str.contains('GL754'), 'object']     = 'Gaia DR3 6663600360560306560'
df.loc[df['object'].str.upper().str.contains('TOI-500'), 'object']     = 'HIP 34269'
df.loc[df['object'].str.upper().str.contains('9CET;BE-CET;9 CET;BECET;HD1835'), 'object']     = '9 CET'
df.loc[df['object'].str.upper().str.contains('HD125248;HD-125248'), 'object']     = 'HD 125248'
df.loc[df['object'].str.upper().str.contains('HD60532;HD060532;HR_2906'), 'object']     = 'HD 60532'
df.loc[df['object'].str.upper().str.contains('TOI-1203'), 'object']     = 'HIP 54779'
df.loc[df['object'].str.upper().str.contains('GL436'), 'object']     = 'HIP 57087'
df.loc[df['object'].str.upper().str.contains('GL536;HD122303'), 'object']     = 'HD 122303'
df.loc[df['object'].str.upper().str.contains('GJ654'), 'object']     = 'HIP 83599'
df.loc[df['object'].str.upper().str.contains('HD133880;HD-133880'), 'object']     = 'HD133880'
df.loc[df['object'].str.upper().str.contains('$EPSILON$,CMA;HD-52089;HD-44743;'), 'object']     = 'EPS CMA'
df.loc[df['object'].str.upper().str.contains('GL628;BD-12 4523'), 'object']     = 'HIP 80824'
df.loc[df['object'].str.upper().str.contains('L-CAR;L CAR'), 'object']     = 'HIP 50847'
df.loc[df['object'].str.upper().str.contains('GL832;GJ 832;GJ-832'), 'object']     = 'HIP 106440'
df.loc[df['object'].str.upper().str.contains('HD-171488;HD171488;V V889 HER;V889-HER'), 'object']     = 'HD 171488'
df.loc[df['object'].str.upper().str.contains('HD211415;HR8501'), 'object']     = 'HD 211415'
df.loc[df['object'].str.upper().str.contains('SW1019-0948;WASP43'), 'object']     = 'Gaia DR3 3767805209112436736'
df.loc[df['object'].str.upper().str.contains('GJ1061;GJ-1061'), 'object']     = 'Gaia DR3 4848140361962951552'
df.loc[df['object'].str.upper().str.contains('GL693'), 'object']     = 'HIP 86990'
df.loc[df['object'].str.upper().str.contains('GJ3135'), 'object']     = 'HIP 9786'
df.loc[df['object'].str.upper().str.contains('TOI-431;HIP26013'), 'object']     = 'HIP 26013'
df.loc[df['object'].str.upper().str.contains('BETA VIR;BETA-VIR;HD-102870;LHS2465_F9V'), 'object']     = 'BET VIR'
df.loc[df['object'].str.upper().str.contains('LRA01_E2_0165'), 'object']     = 'TYC 4799-1733-1'
df.loc[df['object'].str.upper().str.contains('HIP105184;HD202628'), 'object']     = 'HD 202628'
df.loc[df['object'].str.upper().str.contains('49-CNC;49_CNC;HD74521'), 'object']     = '49 CNC'
df.loc[df['object'].str.upper().str.contains('HR9087;HD224926'), 'object']     = 'HD 224926'
df.loc[df['object'].str.upper().str.contains('HD59967;HD  59967;HIP36515;HD059967;HD 59967'), 'object']     = 'HD 59967'
df.loc[df['object'].str.upper().str.contains('GL393;HD851317'), 'object']     = 'HD 851317'
df.loc[df['object'].str.upper().str.contains('55-CANCRI;55CNC;;'), 'object']     = '55 CNC'
df.loc[df['object'].str.upper().str.contains('GJ-479;GL479;GJ479;;'), 'object']     = 'HIP 61629'
df.loc[df['object'].str.upper().str.contains('HR4963;HR4963SPHOT'), 'object']     = 'HR4963'
df.loc[df['object'].str.upper().str.contains('HD204941;HD 204941'), 'object']     = 'HD 204941'
df.loc[df['object'].str.upper().str.contains('WASP-18;WASP18A'), 'object']     = 'HD 10069'
df.loc[df['object'].str.upper().str.contains('GL87;BD+02   348'), 'object']     = 'HIP 10279'
df.loc[df['object'].str.upper().str.contains('HIP-89829;HD 168210'), 'object']     = 'HD 168210'
df.loc[df['object'].str.upper().str.contains('GJ1214'), 'object']     = 'Gaia DR3 4393265392168829056'
df.loc[df['object'].str.upper().str.contains('HD218396;HD218296;V342_PEG'), 'object']     = 'HD 218396'
df.loc[df['object'].str.upper().str.contains('SW1042-0350;SW1042-0350B'), 'object']     = 'Gaia DR3 3778075717162985600'
df.loc[df['object'].str.upper().str.contains('TAU-BOO;TAUBOO'), 'object']     = 'TAU BOO'
df.loc[df['object'].str.upper().str.contains('L98-59;TOI-175;TIO175'), 'object']     = 'TYC 9193-2365-1'
df.loc[df['object'].str.upper().str.contains('HIP-77052;HD140538;PSI SER;HIP77052'), 'object']     = 'HD 140538'
df.loc[df['object'].str.upper().str.contains('HD26967;HD-26967'), 'object']     = 'HD 26967'
df.loc[df['object'].str.upper().str.contains('HD91889;HD091889'), 'object']     = 'HD 91889'
df.loc[df['object'].str.upper().str.contains('DELTA_VEL;DEL-VEL;'), 'object']     = 'DEL VEL'
df.loc[df['object'].str.upper().str.contains('BET DOR;BETA-DOR;BET-DOR'), 'object']     = 'BET DOR'
df.loc[df['object'].str.upper().str.contains('HD77338;HD077338;HD 77338'), 'object']     = 'HD 77338'
df.loc[df['object'].str.upper().str.contains('GJ382;GL382;HIP49986'), 'object']     = 'HIP 49986'
df.loc[df['object'].str.upper().str.contains('GAMMA-SER;'), 'object']     = 'GAM SER'
df.loc[df['object'].str.upper().str.contains('IX_VEL'), 'object']     = 'HIP 40430'
df.loc[df['object'].str.upper().str.contains('GL514'), 'object']     = 'HIP 65859'
df.loc[df['object'].str.upper().str.contains('HIP31592;HD-47205;HIP-31592;HD47205'), 'object']     = 'HD 47205'
df.loc[df['object'].str.upper().str.contains('HD 38677;HD38677'), 'object']     = 'HD 38677'
df.loc[df['object'].str.upper().str.contains('HR1544;HD30739'), 'object']     = 'HD 30739'
df.loc[df['object'].str.upper().str.contains('TOI-1233'), 'object']     = 'HD 108236'
df.loc[df['object'].str.upper().str.contains('HD-179949;HD179949'), 'object']     = 'HD 179949'
df.loc[df['object'].str.upper().str.contains('LP816-60'), 'object']     = 'TYC 6348-400-1'
df.loc[df['object'].str.upper().str.contains('LHS1723'), 'object']     = 'Gaia DR3 3187115498866675456'
df.loc[df['object'].str.upper().str.contains('HD206893;HD 206893'), 'object']     = 'HD 206893'
df.loc[df['object'].str.upper().str.contains('HR5501;HD147513'), 'object']     = 'HD 147513'
df.loc[df['object'].str.upper().str.contains('GL205;GJ205;HD36395'), 'object']     = 'HD 36395'
df.loc[df['object'].str.upper().str.contains('GL752A'), 'object']     = 'HD 173739'
df.loc[df['object'].str.upper().str.contains('EPIC206011496'), 'object']     = 'TYC 5818-486-1'
df.loc[df['object'].str.upper().str.contains('TOI-4399'), 'object']     = 'HIP 94235'
df.loc[df['object'].str.upper().str.contains('HD221420;HD 221420'), 'object']     = 'HD 221420'
df.loc[df['object'].str.upper().str.contains('GL678.1A;'), 'object']     = 'HIP 85665'
df.loc[df['object'].str.upper().str.contains('HD126525;HIP70695;HD 126525'), 'object']     = 'HD 126525'
df.loc[df['object'].str.upper().str.contains('HR1996;HR1996_LIN;HR1996_CIRC;HD38666'), 'object']     = 'HD 38666'
df.loc[df['object'].str.upper().str.contains('SW0415-2206'), 'object']     = 'TYC 5889-271-1'
df.loc[df['object'].str.upper().str.contains('HD117618;HD-117618'), 'object']     = 'HD 117618'
df.loc[df['object'].str.upper().str.contains('HD92499;HD-92499'), 'object']     = 'HD 92499'
df.loc[df['object'].str.upper().str.contains('HD-81009;HR-3724'), 'object']     = 'HD 81009'
df.loc[df['object'].str.upper().str.contains('GJ317;GJ-317;L 675-81'), 'object']     = 'Gaia DR3 5701750715316821248'
df.loc[df['object'].str.upper().str.contains('OMILUP;HD-130807'), 'object']     = 'OMI LUP'
df.loc[df['object'].str.upper().str.contains('C15-0734'), 'object']     = 'HD 137496'
df.loc[df['object'].str.upper().str.contains('HD147513;HD-147513;HIP80337;HD 147513'), 'object']     = 'HD 147513'
df.loc[df['object'].str.upper().str.contains('GL880'), 'object']     = 'HIP 113296'
df.loc[df['object'].str.upper().str.contains('SW0710-3906'), 'object']     = 'TYC 7630-352-1'
df.loc[df['object'].str.upper().str.contains('GJ4100'), 'object']     = 'Gaia DR3 4267422193264246784'
df.loc[df['object'].str.upper().str.contains('GJ3341'), 'object']     = 'Gaia DR3 4827532078086044416'
df.loc[df['object'].str.upper().str.contains('HIP-22263;HIP22263;'), 'object']     = 'HIP 22263'
df.loc[df['object'].str.upper().str.contains('I,CAR;I-CAR'), 'object']     = 'HIP 45101'
df.loc[df['object'].str.upper().str.contains('HR7596;HR-7596'), 'object']     = 'HR 7596'
df.loc[df['object'].str.upper().str.contains('HD117207;HD 117207'), 'object']     = 'HD 117207'
df.loc[df['object'].str.upper().str.contains('BD+20  2465;GL388;AD-LEO;HD901741'), 'object']     = 'TYC 1423-174-1'
df.loc[df['object'].str.upper().str.contains('HD-188041;HD188041'), 'object']     = 'HD 188041'
df.loc[df['object'].str.upper().str.contains('GL410;HD95650;HD095650'), 'object']     = 'HD 95650'
df.loc[df['object'].str.upper().str.contains('HD67200;HD 67200'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('GL358'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD52265;HD  52265'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('GJ1132'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD76700;HD-76700'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('ALP-SCL;HD5737'), 'object']     = 'ALF SCL'
df.loc[df['object'].str.upper().str.contains('V V1358 ORI;HD-43989;HD43989;HD043989'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('GL213'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD 203384;HD203384'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('GJ3470;LP 424-4;HD85512'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('EPIC249622103'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD-116114'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD204313;HD 204313'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('TOI-125;TOI_125'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('TOI-2494;TIC282576340'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('DSTUC;HD 222259;HD222259'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD28471;HD 28471'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('ALF-GRU'), 'object']     = 'ALF GRU'
df.loc[df['object'].str.upper().str.contains('HIP-54582;HD97037;HIP54582;HD 97037'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD-128429'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD013445_GJ86_HR637;HD 13445;HD013445'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD89839;HD 89839'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('BD-061339;GJ221'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD 219077;HD219077'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD-201601;GAMMA-EQU;HD201601'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('EPIC228801451'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('TOI-1054'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HIP102040;HD 197076'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD 92788;HD92788'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD-27536;HD27536'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD110014;HD-110014;CHI-VIR_K2III'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD216437;HD 216437'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD34816;HD_34816_4.29'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD114729;HD 114729'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('GJ-729;HD000011;GL729'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD42936;HD 42936'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD114783;HD 114783'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HR-1654;HD-32887'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HIP-30476;HD45289;HIP30476'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD099211;HD99211'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HN-PEG;HD 206860;HD206860'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD43162;HD  43162'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('BD-082823'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD029391;HD29391'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD564;HD 564'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD-3405;HD3405'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HIP-11915;HD16008;HIP11915'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('GJ2066'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('TOI-544;TIC 713009339'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('BD-210397'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD4308;HD221356'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD55;HD     55'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('TOI179;HD018599;TOI-179;TIO179'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('TOI-421'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('GJ361'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD218994;HD218994A'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('AI_VEL'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD215641;HD 215641'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('SW1233-1008;SW1233-1003'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('TOI-1130'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD33142;HD033142'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('TOI-220'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HIP-29432;HIP29432;HD42618'), 'object']     = 'HD 42618'
df.loc[df['object'].str.upper().str.contains('EPS-CEN;HD-118716'), 'object']     = 'EPS CEN'
df.loc[df['object'].str.upper().str.contains('GAM-VEL;HD-68243'), 'object']     = 'GAM VEL'
df.loc[df['object'].str.upper().str.contains('2M1130+0735'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('2M1129-0127'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('SW2359-3501'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('SW0604-1658'), 'object']     = ''
df.loc[df['object'].str.upper().str.contains('HD128620;HD128620A;HR5459;ALPHACENA;ALPHA-CEN-A;ALPCENA'), 'object']     = 'ALF CEN A'
df.loc[df['object'].str.upper().str.contains('HD128621;ALF CEN B;HR5460;HR-5460;ALPHACENB;ALPCENB'), 'object']     = 'ALF CEN B'
df.loc[df['object'].str.upper().str.contains('ALPHA-CEN;ALPHA-CENTAURI'), 'object']     = 'ALF CEN'
#df.loc[df['object'].str.upper().str.contains(''), 'object']     = ''


mask_sun=df['object'].str.lower()=='sun'
mask_ssbodies=df['object'].str.lower().str.contains('moon|ceres|vesta|astraea|hebe|europa|ganymede|herculina|irene|iris|juno|jupiter|massalia|melpomene|parthenope|pallas|venus|73p-c|titan') | (df['object'].str.lower() == 'io')
mask_test=(df['ra']<0) | (df['ra']>360) | (df['dec']<-90) | (df['dec']>30) | df['object'].str.lower().str.contains('zenit|preset|test|focus|stab')

mask_solsys=mask_sun|mask_ssbodies

sun        = df[mask_sun]
ssbodies   = df[mask_ssbodies]
solsys     = df[mask_solsys]
remaining  = df[~mask_solsys]
testdata   = remaining[mask_test]
astrodata  = remaining[~mask_test]

#sun.to_csv('output/sun.csv', index=False)
#ssbodies.to_csv('output/solssys_bodies.csv', index=False)
#testdata.to_csv('output/test_data.csv', index=False)
#astrodata.to_csv('output/astro_data.csv', index=False)


mask_clusters= (astrodata['prog_title'].str.lower().str.contains('cluster')|astrodata['prog_title'].str.lower().str.contains('m67')|astrodata['object'].str.lower().str.startswith('ngc') |astrodata['object'].str.lower().str.startswith('ic')|astrodata['object'].str.lower().str.startswith('47tuc')|astrodata['object'].str.lower().str.startswith('m67')|astrodata['object'].str.lower().str.startswith('melotte')|astrodata['object'].str.lower().str.contains('terzan'))
mask_corot   = astrodata['object'].str.lower().str.startswith('ira0') | astrodata['object'].str.lower().str.startswith('sra0') | astrodata['object'].str.lower().str.startswith('src0') | astrodata['object'].str.lower().str.startswith('lra0') | astrodata['object'].str.lower().str.startswith('lrc0')
mask_kepler2 = astrodata['object'].str.lower().str.startswith('k2') | astrodata['object'].str.lower().str.startswith('epic')
mask_tess    = astrodata['object'].str.lower().str.startswith('toi') | astrodata['object'].str.lower().str.startswith('tic')
mask_hat     = astrodata['object'].str.lower().str.startswith('hats') | astrodata['object'].str.lower().str.startswith('hat')
mask_swasp   = astrodata['object'].str.lower().str.startswith('sw') | astrodata['object'].str.lower().str.contains('wasp')
mask_ogle    = astrodata['object'].str.lower().str.startswith('ogle') | astrodata['object'].str.lower().str.startswith('blg')| astrodata['object'].str.lower().str.startswith('lmc')| astrodata['object'].str.lower().str.startswith('smc')
mask_all     = mask_tess | mask_clusters | mask_corot | mask_hat | mask_swasp | mask_ogle | mask_kepler2


tess     =  astrodata[mask_tess]
clusters =  astrodata[mask_clusters]
corot    =  astrodata[mask_corot]
kepler2  =  astrodata[mask_kepler2]
hat      =  astrodata[mask_hat]
swasp    =  astrodata[mask_swasp]
ogle     =  astrodata[mask_ogle]
stelle   =  astrodata[~mask_all]

clusters.to_csv('output/clusters.csv', index=False)
corot.to_csv('output/corot.csv', index=False)
kepler2.to_csv('output/kepler2.csv', index=False)
tess.to_csv('output/tess.csv', index=False)
hat.to_csv('output/hat.csv', index=False)
swasp.to_csv('output/swasp.csv', index=False)
ogle.to_csv('output/ogle.csv', index=False)
stelle.to_csv('output/stelle.csv', index=False)


#clustering(astrodata,'astrodata',0.01)

sys.exit()



nums={'all obs':len(df),'sun':len(sun),'solsys':len(ssbodies),'test':len(testdata),'astro':len(astrodata),'clusters':len(clusters),'corot':len(corot),'kepler2':len(kepler2),'tess':len(tess),'hat':len(hat),'swasp':len(swasp),'ogle':len(ogle),'stars':len(stelle)}

for key, value in nums.items():
    print(f"{key}: {value}")


dataframes={'stelle':stelle,'tess':tess,'clusters':clusters,'corot':corot,'kepler2':kepler2,'hat':hat,'swasp':swasp,'ogle':ogle,'ssbodies':ssbodies,'sun':sun,'testdata':testdata}

for key, df in dataframes.items():
    df['OBSGROUP'] = key

#for name, df in dataframes.items():
#    #print(name)
#    clustering(df,name,0.01)


