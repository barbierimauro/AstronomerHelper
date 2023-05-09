SELECT distinct instrument,min(period),max(period),min(date_obs),max(date_obs), sum(exposure)/86400  as ttot, --sum(exposure)/31536000  as ttoty 
FROM dbo.raw
where dp_cat like 'science%'
group by instrument
order by ttot desc
