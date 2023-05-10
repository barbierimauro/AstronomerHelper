#!/bin/bash
# 2022 Mauro Barbieri
# remove exif tags from photos and rescale to 50% of original size and convert to png

ext=JPG

ls *.$ext | while read f
do
	opng="${f/.$ext/.png}"
        echo $opng
        exiftool -all= $f -o tmp.jpg
        convert -resize 50% tmp.jpg $opng
        rm -rf tmp.jpg
done
