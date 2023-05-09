#!/bin/bash
#2022 Mauro Barbieri
#make a mp4 movie form a sequence of images

# the script assumes that all the files have the same shape (MANDATORY)
# and that their name are in a alphanumerical order

# no check are performed on the size of the images
# any previous output is deleted


conv=true # create the images with sequential names
nrep=10 # how many times the last frame is repeated

#1 
# convert all the files in the actual directory in to a sequence 
# of files with the proper names for ffmpeg

if [ $conv = true ]; then
    rm imgtmp__*.jpg
    for p in `ls *jpg`
    do
	let i=$i+1
	n=$(printf "%04d" $i)
        o=imgtmp__$n.jpg
	ln -s $p $o
        # make some geometric operation with convert
        # in this case the previous operation with ln need to be substituted by cp
	#convert $p -crop +0+1 $0 
    done 
    if [ $nrep -gt 0 ]; then
	o1=$o
    	seq 1 $nrep | while  read j
    	do
		let i=$i+1
		n=$(printf "%04d" $i)
        	o=imgtmp__$n.jpg
		ln -s $o1 $o
    	done
    fi
fi

#2
# make the movie 

rm -rf out.mp4
ffmpeg -framerate 2 -i imgtmp__%04d.jpg -c:v libx264 -pix_fmt yuv420p out.mp4

#rm -rf imgtmp__*.jpg

