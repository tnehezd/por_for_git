#!/bin/bash

out=kep
index=0

#for ((i=1; i<=9999; i++))
#   do
   cat kep.list | \
   while read kep dummy ; do
        if test -f $kep.dat ; then
#		echo "van"
                gnuplot -e "set terminal jpeg; set xrange[-15:15]; set yrange[-15:15]; plot '$kep.dat' u 2:3" > $out$index.jpeg
		let index+=1
	else 
		echo $kep_2.dat "nincs"
	fi
    done
	
#done

ffmpeg -f image2 -i kep%d.jpeg movie2_2.mpeg

