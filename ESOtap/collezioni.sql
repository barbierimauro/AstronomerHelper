SELECT s_ra as ra, s_dec as dec, snr, (t_max+t_min)/2 as mjd, t_exptime as exptime, em_xel as npt, em_res_power as respow, em_min as wlmin, em_max as wlmax, target_name as target, o_calib_status as calib, dp_id
FROM ivoa.obscore
where 
obs_collection like 'giraffe%'
and snr>3
