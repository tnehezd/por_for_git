#!/bin/bash

TIME=time.dat

NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "100000. 250 $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.6.c -o $INIT -lm
./$INIT  -ri 1. -ro 40. -rdzei 2.7 -drdzei 2. -rdzeo 20. -drdzeo 2. -index 0.5 -alpha 0.01 -n 500 -onesize 0.1 -md 0.02

DRIFT=drift$NOW

gcc -O3 dustdrift1Dver16_rdf.c -o $DRIFT -lm
./$DRIFT -tStep 0.5 -n 500 -twopop 0 
wait

