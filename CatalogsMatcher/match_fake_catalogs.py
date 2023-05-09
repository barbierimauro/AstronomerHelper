#!/usr/bin/python3
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

def calculate_metric(phys_ref, phys_cmp, err_phys_ref):
    if pd.notnull(phys_ref) and pd.notnull(phys_cmp) and abs(err_phys_ref) > 0:
        metric_n = 1
        metric = ((phys_ref - phys_cmp) / err_phys_ref) ** 2
    else:
        metric_n = 0
        metric = 0
    return metric, metric_n

def cross_match(ref_catalog, cmp_catalog):
    combined_columns = ref_catalog.columns.tolist() + cmp_catalog.columns.tolist() + ['metric', 'metric_plx', 'metric_pmra', 'metric_pmde', 'metric_rv', 'metric_vmag', 'metric_vicol', 'angular_separation']
    matches = pd.DataFrame(columns=combined_columns)

    for index, ref_star in ref_catalog.iterrows():
        print(index)
        ref_coord = SkyCoord(ra=ref_star['ra']*u.degree, dec=ref_star['de']*u.degree)
        e_pos_ref = np.sqrt(ref_star['e_ra']*ref_star['e_de'])/1000.0 # in arcsecond

        for _, cmp_star in cmp_catalog.iterrows():
            cmp_coord = SkyCoord(ra=cmp_star['ra']*u.degree, dec=cmp_star['de']*u.degree)
            angular_separation = ref_coord.separation(cmp_coord).arcsecond

            if angular_separation <= 10:
                dist = 0
                num = 0
                metric_values = []

                # Calculate metric_0 for angular separation
                metric_0 = angular_separation / e_pos_ref
                dist += metric_0
                num += 1
                metric_values.append(metric_0)

                for physical_var in ['plx', 'pmra', 'pmde', 'rv', 'vmag', 'vicol']:
                    metric, metric_n = calculate_metric(ref_star[physical_var], cmp_star[physical_var], ref_star['e_' + physical_var])
                    dist += metric
                    num += metric_n
                    metric_values.append(metric)

                if num > 0:
                    dist /= num

                match_data = pd.concat([ref_star, cmp_star])
                match_data['metric'] = dist
                match_data['metric_plx'] = metric_values[1]
                match_data['metric_pmra'] = metric_values[2]
                match_data['metric_pmde'] = metric_values[3]
                match_data['metric_rv'] = metric_values[4]
                match_data['metric_vmag'] = metric_values[5]
                match_data['metric_vicol'] = metric_values[6]
                match_data['angular_separation'] = angular_separation
                matches = matches.append(match_data, ignore_index=True)
                #print(match_data)
                

    return matches


# Example usage
ref_catalog = pd.read_csv('catalog_ref.csv')
cmp_catalog = pd.read_csv('catalog_cmp.csv')

matches = cross_match(ref_catalog, cmp_catalog)

# Save DataFrame to CSV and FITS files
matches.to_csv('matched_catalog.csv', index=False)
#astropy_table = Table.from_pandas(output_df)
#astropy_table.write('matched_catalog.fits', overwrite=True)
