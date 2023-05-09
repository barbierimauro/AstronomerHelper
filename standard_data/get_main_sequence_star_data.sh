#!/bin/bash
msfile="EEM_dwarf_UBVIJHK_colors_Teff.txt"
url="https://www.pas.rochester.edu/~emamajek/"$msfile

# Download the file if it doesn't exist
if [ ! -e $msfile ]; then
   #echo "File does not exist, downloading it now..."
   wget $url
else
   true
   #echo "File exists, no need to download."
fi

out1="msdata.tmp"
touch $out1

readmefile="main_sequence_data.readme"
touch $readmefile


awk '
BEGIN {
  spt_found = 0
}
{
  if ($0 ~ /#SpT/) {
	if (spt_found == 0) {
  	print > "'"$out1"'"
  	spt_found = 1
	} else {
  	spt_found = 0
	}
  } else if (spt_found == 1) {
	print > "'"$out1"'"
  } else {
	print > "'"$readmefile"'"
  }
}' "$msfile"


sed -i '/^#SpT/ { s/R_Rsun/Rsun/g; s/-//g; s/Ks/K/g; s/Ic/I/g; s/:/ /g; s/  */,/g; s/^,//g; s/,$//g; }' $out1
sed -e 's/\:/ /g' -e 's/\.\.\.\.\./\.\.\.  /g' -e 's/\.\.\.\./\.\.\. /g' -e 's/Msun,#SpT/Msun/g' -e 's/#//g' $out1 > _oo
mv _oo $out1

chmod +x get_main_sequence_star_data_fortran_reformat
chmod +x get_main_sequence_star_data_python_converter.py
./get_main_sequence_star_data_fortran_reformat
./get_main_sequence_star_data_python_converter.py
chmod -x get_main_sequence_star_data_fortran_reformat
chmod -x get_main_sequence_star_data_python_converter.py

rm -rf msdata.tmp msdata.csv

echo
echo "README:  main_sequence.readme"
echo "DATA  :  main_sequence_data.csv"
echo
