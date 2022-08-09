#!/bin/bash
inp=$1
pref="${inp/.csv/}"
odir=$pref
nrow=$2
cwd=$(pwd -P)

narg=$#
if [ $narg -eq 2 ]; then
   true
else
   echo "ERROR the number of argument must be exactly 2:"
   echo "      the first argument is the csv file name"
   echo "      the second argument is the number of lines in which you want to split"
   echo
   echo "syntax:  splitcsv dataset.csv 10000"
   exit
fi

ifile=$cwd"/"$inp
rpath=$(readlink -f "$inp")
if [ $ifile = $rpath ]; then
   #echo uguali
   true
else
   echo "ERROR the input file is expect to be in the current working directory"
   echo "input file " $inp
   echo "realpath   " $rpath
   exit
fi


kinp=$(file $inp | grep -c "CSV")
kint=$(file $inp)
if [ $kinp -eq 1 ]; then
   true
else
   echo "ERROR the input file is not reported as a CSV file (or doesn't exist)"
   echo $kint
   exit
fi


if [ $nrow -eq $nrow ] 2> /dev/null
then
   if [ $nrow -gt 1 ]; then
      true
      if [ $nrow -lt 100 ]; then
         echo "WARNING: the number of lines is below 100"
      fi
   else
      echo "ERROR the second argument need to be an integer larger than 1"
      echo "nrow="$nrow
      exit
   fi
else
   echo "ERROR the second argument need to be an integer larger than 1"
   echo "nrow="$nrow
   exit
fi





head -1 $inp > hdr
hdr=$cwd"/hdr"



rm -rf $odir
mkdir $odir
cd $odir
echo splitting $inp
fwfl=filename_without_first_line
cat $cwd"/"$inp | sed 1d > $fwfl
split -u -e -l $nrow --numeric-suffixes=1 -a 3 --additional-suffix=".tmpcsv" $fwfl $pref"."

i=0
nf=$(ls *.tmpcsv|wc -l)
for f in *.tmpcsv
do
  let i=$i+1
  if [ $nrow -gt 100000 ]; then
     echo "making slice n. "$i "of "$nf
  fi
  ofile="${f/.tmpcsv/.csv}"
  cat $hdr $f >$cwd"/"$ofile
  rm -rf $f 2> /dev/null
done

cd $cwd
rm -rf $hdr $odir 2> /dev/null

echo "splitted files are availables in "$cwd
echo


