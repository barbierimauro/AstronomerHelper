#!/bin/bash
#2022 Mauro Barbieri
#make a mp4 movie form a sequence of pair of images

# the script take two images with same size, and sticht together horizontally or vertically
# then rename the files in a format compatible with ffmpeg and make a movie

# the script assumes that all the files in each directory have the same shape (MANDATORY)
# and that their name are in a alphanumerical order

# no check are performed on the size and the number of the images
# any previous output is deleted


conv=true # create the images with sequential names
nrep=10 # how many times the last frame is repeated
dir1=raw
dir2=red

#rm -rf tmp
#mkdir tmp
#mkdir ovideo

if [ $conv = true ]; then
    rm -rf ovideo/appendimg*.jpg
    cd $dir1
    for p in `ls *jpg`
    do
	echo $p
	let i=$i+1
	n=$(printf "%04d" $i)
	cp $p ../tmp/1.jpg
	cp ../$dir2/$p ../tmp/2.jpg
	cd ../tmp
        pwd
	convert -quality 100 +append *jpg out.jpg # stitch horizontally
	#convert -quality 100 -append *jpg out.jpg # stitch vertically
	rm -rf 1.jpg 2.jpg
        o="appendimg"$n".jpg"
	mv out.jpg ../ovideo/$o
	cd ../$dir1
    done
    cd ..
    if [ $nrep -gt 0 ]; then
	o1=$o
    	seq 1 $nrep | while  read j
    	do
		let i=$i+1
		n=$(printf "%04d" $i)
        	o="imgtmp__"$n".jpg"
		cp ovideo/$o1 ovideo/$o
    	done
    fi
fi


cd ovideo
ffmpeg -framerate 2 -i appendimg%04d.jpg -c:v libx264 -pix_fmt yuv420p out.mp4

#rm -rf tmp


