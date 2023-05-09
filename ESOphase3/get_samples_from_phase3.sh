#!/bin/bash
# Download some "random" ADP files from phase3 archive
# It download nf files per each combination of dataproduct_type and dataproduct_subtype

author="Mauro Barbieri"
year=2022
email=maurobarbieri.science@gmail.com


#combinations of dataproduct_type and dataproduct_subtype
#image,null
#spectrum,null
#visibility,null
#cube,ifs
#image,exposure
#image,fluxmap
#image,pawprint
#image,tile
#measurements,catalog
#measurements,catalogtile
#measurements,srctbl
#image,deeppawprint
#image,deeptile

#example TAP query
#SELECT top 1 dp_id FROM ivoa.ObsCore where dataproduct_type like 'image' and dataproduct_subtype like 'deep pawprint' and p3orig is not null and dp_id is not null


nf=3

url01="http://archive.eso.org/tap_obs/sync?REQUEST=doQuery&LANG=ADQL&MAXREC=100&FORMAT=csv&QUERY=SELECT%20top%20"$nf"%20dp_id%20FROM%20ivoa.ObsCore%20where%20dataproduct_type%20like%20%27"
urla2="%27%20and%20dataproduct_subtype%20is%20null%20and%20p3orig%20is%20not%20null%20and%20dp_id%20is%20not%20null"
urlb2="%27%20and%20dataproduct_subtype%20like%20%27"
urlb3="%27%20and%20p3orig%20is%20not%20null%20and%20dp_id%20is%20not%20null"


ark="https://dataportal.eso.org/dataportal_new/file/"

mkdir p3samples 2>/dev/null
cd p3samples
i=0
awk -F, '{printf"%s %s\n",$1,$2}' ../list_dataproducts.csv | while read d1 d2
do
	let i=$i+1
    mkdir -p $d1/$d2 2>/dev/null
	cd $d1/$d2
	if (( $i >= 1     )) && (( $i <= 3 )); then
	   qsql=$url01$d1$urla2
	elif (( $i > 3    )) && (( $i <= 11 )); then
	   qsql=$url01$d1$urlb2$d2$urlb3
	elif (( $i ==  12 )); then
	   #"deep pawprint"
	   qsql="http://archive.eso.org/tap_obs/sync?REQUEST=doQuery&LANG=ADQL&MAXREC=100&FORMAT=csv&QUERY=SELECT%20top%20"$nf"%20dp_id%20FROM%20ivoa.ObsCore%20where%20dataproduct_type%20like%20%27image%27%20and%20dataproduct_subtype%20like%20%27deep%20pawprint%27%20and%20p3orig%20is%20not%20null%20and%20dp_id%20is%20not%20null"
	elif (( $i ==  13 )); then
	   #"deep tile"
	   qsql="http://archive.eso.org/tap_obs/sync?REQUEST=doQuery&LANG=ADQL&MAXREC=100&FORMAT=csv&QUERY=SELECT%20top%20"$nf"%20dp_id%20FROM%20ivoa.ObsCore%20where%20dataproduct_type%20like%20%27image%27%20and%20dataproduct_subtype%20like%20%27deep%20tile%27%20and%20p3orig%20is%20not%20null%20and%20dp_id%20is%20not%20null"
	fi
    wget -O tmp__ $qsql 2>/dev/null
    grep -v dp_id tmp__ | while read dpid
    do
		echo $i $d1 $d2 $dpid
		dpidf=$dpid".fits"
		wget -O $dpid "$ark""$dpid" 2>/dev/null
		mv $dpid $dpidf
		##file $dpidf
    done
    rm -rf tmp__
    cd ../..
    echo
done

exit







