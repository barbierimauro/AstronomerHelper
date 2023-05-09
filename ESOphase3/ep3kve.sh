#!/bin/bash
# ep3kve : eso phase 3 fits keyword values extractor
#
author="Mauro Barbieri"
year=2022
email=maurobarbieri.science@gmail.com
ver=2022-10-06T10:30:55_v0.1
vern=${ver:21:23}
verd=${ver:0:19}

function logo() {
echo "+---------------------------------+"
echo "| ep3kve                          |"
echo "| ver:" $vern"                        |"
echo "|                                 |"
echo "| "$author "(c)" $year   "        |"
echo "| "$email                        "|"
echo "+---------------------------------+"
#echo "This is free software; see the source for copying conditions.  There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE."
echo
}

logo 

kwlist=(SIMPLE DATE BITPIX NAXIS NAXIS1 NAXIS2 NAXIS3 EXTEND XTENSION EXTNAME PRODCATG ASSOC1 ASSOC2 ASSOC3 ASSON1 ASSON2 ASSON3 ASSOM1 ASSOM2 ASSOM3 ORIGIN TELESCOP TELESC1 TELESC2 TELESC3 INSTRUME INSTR1 INSTR2 INSTR3 FILTER FILTER1 FILTER2 FILTER3 OBJECT RA DEC EQUINOX RADESYS RADECSYS TIMESYS EXPTIME TEXPTIME MJD-OBS MJD-END PROG_ID PROGID1 PROGID2 PROGID3 OBID1 OBID2 OBID3 OBID4 OBID5 OBID6 OBID7 OBID8 OBID9 NCOMBINE OBSTECH FLUXCAL PROCSOFT REFERENC PROV1 PROV2 PROV3 PROV4 PROV5 PROV6 PROV7 PROV8 PROV9 PROVXTN BUNIT GAIN DETRON EFFRON WEIGHT CRVAL1 CRVAL2 CRVAL3 CRPIX1 CRPIX2 CRPIX3 CTYPE1 CTYPE2 CTYPE3 CUNIT1 CUNIT2 CUNIT3 CD1_1 CD1_2 CD1_3 CD2_1 CD2_2 CD2_3 CD3_1 CD3_2 CD3_3 CSYER1 CSYER2 CSYER3 CRDER1 CRDER2 CRDER3 PHOTZP PHOTZPER PHOTSYS SPECSYS EXT_OBJ CONTNORM TOT_FLUX FLUXERR WAVELMIN WAVELMAX LAMRMS LAMNLIN SPEC_BIN SPEC_ERR SPEC_SYE RA_ERR DEC_ERR NELEM VOCLASS VOPUB TITLE APERTURE TELAPSE TMID SPEC_VAL SPEC_BW BNOISE MAPMODE FEBE1 FEBE2 FEBE3 CONTENT INSMODE BASE_MIN BASE_MAX NUM_CHAN VIS2ERR T3PHIERR STOKES HDUCLASS HDUCLAS1 HDUCLAS2 HDUCLAS3 HDUCLAS4 HDUCLAS5 HDUCLAS6 HDUCLAS7 HDUCLAS8 HDUCLAS9 HDUCLAS10 HDUDOC HDUVERS SCIDATA ERRDATA QUALDATA CONFDATA BKGDATA BKGERR BKGCONF ABMAGLIM PIXNOISE MAGLIM1 MAGLIM2 MAGLIM3 MAGLIM4 MAGLIM5 MAGLIM6 MAGLIM7 MAGLIM8 MAGLIM9 ABMAGSAT PSF_FWHM ELLIPTIC SNR SPEC_RES SKY_RES SKY_RERR STREHL ARCFILE CHECKSUM DATASUM ORIGFILE P3ORIG NDIT NJITTER NOFFSETS NUSTEP FPRA1 FPRA2 FPRA3 FPRA4 FPRA5 FPRA6 FPRA7 FPRA8 FPDE1 FPDE2 FPDE3 FPDE4 FPDE5 FPDE6 FPDE7 FPDE8 SKYSQDEG M_EPOCH APMATCHD TXLNK1 TXLNK2 TXLNK3 TXRGF TXCTY NOESODAT TFIELDS TTYPE1 TTYPE2 TTYPE3 TTYPE4 TTYPE5 TTYPE6 TTYPE7 TTYPE8 TTYPE9 TTYPE10 TFORM1 TFORM2 TFORM3 TFORM4 TFORM5 TFORM6 TFORM7 TFORM8 TFORM9 TFORM10 TCOMM1 TCOMM2 TCOMM3 TCOMM4 TCOMM5 TCOMM6 TCOMM7 TCOMM8 TCOMM9 TCOMM10 TUNIT1 TUNIT2 TUNIT3 TUNIT4 TUNIT5 TUNIT6 TUNIT7 TUNIT8 TUNIT9 TUNIT10 TUTYP1 TUTYP2 TUTYP3 TUTYP4 TUTYP5 TUTYP6 TUTYP7 TUTYP8 TUTYP9 TUTYP10 TUCD1  TUCD2  TUCD3  TUCD4  TUCD5  TUCD6  TUCD7  TUCD8  TUCD9  TUCD10  TDMIN1 TDMIN2 TDMIN3 TDMIN4 TDMIN5 TDMIN6 TDMIN7 TDMIN8 TDMIN9 TDMIN10 TDMAX1 TDMAX2 TDMAX3 TDMAX4 TDMAX5 TDMAX6 TDMAX7 TDMAX8 TDMAX9 TDMAX10 TNULL1 TNULL2 TNULL3 TNULL4 TNULL5 TNULL6 TNULL7 TNULL8 TNULL9 TNULL10 TZERO1 TZERO2 TZERO3 TZERO4 TZERO5 TZERO6 TZERO7 TZERO8 TZERO9 TZERO10 TSCAL1 TSCAL2 TSCAL3 TSCAL4 TSCAL5 TSCAL6 TSCAL7 TSCAL8 TSCAL9 TSCAL10 EXTNAME EXTVER EXTLEVEL COMMENT)

debug=0

BDIR_base="/diskb/phase3data/ftp/programs/"


progid="$1"
WDIR="/home/phase3/mbarbier/WORK""${progid}"
BDIR="${BDIR_base}""${progid}"
mkdir -p $WDIR



function exthdr() {
echo processing $BDIR
echo output directory $WDIR
echo
cd $BDIR

nf=$(ls *fits | wc -l)
i=0

ls *fits | while read f
do
   let i=$i+1
   echo $f file $i of $nf
   nend=$(dfits -x 0 $f | grep -c ^END)
   let next=$nend-1
   echo "total extensions = "$next
   echo extracting main header ...
   ohdr=$WDIR"/"$f"_ext_0.hdr"
   dfits $f > $ohdr
   if [ $next -ge 1 ]; then
      seq 1 $next | while read i
      do
        echo extracting header extension $i ...
        ohdr=$WDIR"/"$f"_ext_"$i".hdr"
        dfits -x $i $f > $ohdr
      done
   fi
   echo
done
}


function extkw() {
cd $WDIR
out=report.csv

str="filename"
for kname in "${kwlist[@]}";
do
  str="${str}"",""${kname}"
done
echo $str>$out

i=0
nf=$(ls *hdr | wc -l)


ls *hdr | while read f
do
  str="${f}"
  let i=$i+1
  echo $f file $i of $nf
  for kname in "${kwlist[@]}";
  do
    kname_val=$(cat $f | grep -m1 ^"${kname}"| cut -d= -f2 | cut -d/ -f1 | sed -e "s/'//g" -e "s/ //g")
    if [ $debug -eq 1 ]; then
       mm=$(cat $f | cut -d= -f1 | sed -e "s/ //g" | grep -c ^"${kname}"$)
       nn=$(cat $f | cut -d= -f1 | grep ^"${kname}")
       echo $kname = $kname_val, $mm, $nn
    fi
    str="${str}"",""${kname_val}"
   done
   echo $str>>$out
done
echo
}


exthdr
extkw

stilts  tpipe ifmt=csv in=report.csv cmd='stats Index Name Class Shape Elsize Units Description UCD UCD_desc NGood NBad Cardinality Minimum Maximum MinPos MaxPos Sum Mean Median StDev Variance Skew Kurtosis Quartile1 Quartile2 Quartile3' ofmt=csv out=report.stats.csv

echo done

