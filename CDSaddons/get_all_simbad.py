#!usr/bin/python3
from astroquery.utils.tap.core import TapPlus
from astropy.io import fits

# Create a connection to the TAP service
tap_service = TapPlus(url="http://simbad.u-strasbg.fr/simbad/sim-tap")

# Define the base of your query
query_base = """
SELECT main_id,otype,otypes,sp_type,morph_type,nbref,ra,dec,B,V,J,H,K,plx_value,pmra,pmdec,rvz_radvel,rvz_redshift,rvz_err,avg(teff) as Teff,avg(log_g) as logg,
avg(fe_h) as FeH,galdim_majaxis,galdim_minaxis,galdim_angle,
coo_err_maj,coo_err_min,coo_err_angle,plx_err,coo_qual,coo_wavelength,rvz_type,sp_qual,
oid,hpx,ids,update_date
FROM basic
left outer JOIN allfluxes ON basic.oid=allfluxes.oidref
left outer JOIN alltypes ON basic.oid=alltypes.oidref
left outer JOIN ids on basic.oid=ids.oidref
left outer JOIN mesFe_h on basic.oid=mesFe_h.oidref
WHERE oid >= {} AND oid < {}
GROUP BY main_id,otype,otypes,sp_type,morph_type,nbref,ra,dec,B,V,J,H,K,plx_value,pmra,pmdec,rvz_radvel,rvz_redshift,rvz_err,galdim_majaxis,galdim_minaxis,galdim_angle,coo_err_maj,coo_err_min,coo_err_angle,plx_err,coo_qual,coo_wavelength,rvz_type,sp_qual,oid,hpx,ids,update_date
"""

# Loop over the 20 queries
for i in range(20):
    lower_bound = i * 1e6
    upper_bound = (i + 1) * 1e6

    # Format the query with the current bounds
    query = query_base.format(lower_bound, upper_bound)

    # Execute the query
    result = tap_service.launch_job(query).get_results()

    # Save the results in a FITS file
    result.write(f"results_{i}.fits", format="fits")
