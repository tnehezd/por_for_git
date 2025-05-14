#!/bin/sh

list=`ls *.list`

for i in ${list} ; do 
    echo $i 
    
#    if test 
#    then
	cat $i | grep -i $1 
	cat $i | grep -i $2
#    fi
done

