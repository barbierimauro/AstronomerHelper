--SELECT distinct instrument, filter_path, count(*) as n
--select count(*)
--select object,target,date_obs,mjd_obs,exposure,telescope,instrument,filter_path,ra,dec,gal_lon,gal_lat,ecl_lon,ecl_lat,tel_az,tel_alt,tel_ambi_pres_start,tel_ambi_rhum,tel_ambi_fwhm_start,dp_cat,dp_tech,dp_type,ins_mode,ob_name,tpl_id,origfile,period,prog_id,prog_type,obs_mode,release_date,dp_id
select count(*)
FROM dbo.raw
where (dp_cat like 'science' or dp_cat like 'calib')
and dp_tech like 'ima%'
and dp_tech not like '%flat%'
and dp_tech not like '%dark%'
and dp_tech not like '%bias%'
and dp_tech not like '%dome%'
and dp_tech not like '%lamp%'
and dp_type not like '%flat%'
and dp_type not like '%dark%'
and dp_type not like '%bias%'
and dp_type not like '%dome%'
and dp_type not like '%lamp%'
and exposure > 1
and origfile not like '%spec%'
and origfile not like '%dark%'
and origfile not like '%flat%'
and origfile not like '%bias%'
and (grat_path is null or grat_path like 'free')
and (gris_path is null or gris_path like 'free')
and (slit_path is null or slit_path like 'free')
and (instrument like 'efosc%'
or instrument like 'emmi%'
or instrument like 'fors1'
or instrument like 'fors2'
or instrument like 'grond'
or instrument like 'hawki'
or instrument like 'isaac'
or instrument like 'naco'
or instrument like 'omegacam'
or instrument like 'sofi%'
or instrument like 'susi%'
or instrument like 'timmi%'
or instrument like 'vimos'
or instrument like 'vircam'
or instrument like 'visir'
or instrument like 'wfi'
or instrument like 'xshooter')
and release_date < '2022-07-01'
and ra_pnt between 0 and 360
and dec_pnt between -90 and 45
and ra between 0 and 360
and dec between -90 and 45
and tel_alt between 10 and 90
and tel_az between 0 and 360
--group by instrument, filter_path
--order by n desc
