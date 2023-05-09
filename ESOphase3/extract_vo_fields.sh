#!/bin/bash

if [ -z "$1" ]; then
  echo "Error: no input file provided"
  exit 1
fi

input="$1"

if [ ! -e "$input" ]; then
  echo "Error: file '$input' does not exist"
  exit 1
fi

if [ ! -s "$input" ]; then
  echo "Error: file '$input' is empty"
  exit 1
fi

if ! file "$input" | grep -q "FITS"; then
  echo "Error: file '$input' is not a FITS file"
  exit 1
fi

filename=$(basename "$input" .fits)

echo "preparing $input"
dfits -x 0 $input > "$filename.hdr"
grep -E "TTYPE|TFORM|TCOMM|TUNIT|TUTYP|TUCD|TDMINI|TDMAX|TNULL" "$filename.hdr" > keywords_vo.list
echo
echo "extract VO fields"
./step1.extract_vo_fields.py
echo
echo "check physical units"
./step2.check_vo_units.py
echo
echo "check ucd values"
./step3.check_vo_ucds.py
echo
