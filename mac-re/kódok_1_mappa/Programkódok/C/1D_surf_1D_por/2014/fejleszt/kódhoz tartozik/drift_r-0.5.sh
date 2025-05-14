#!/bin/bash

TIME=time.dat

# ALPHA 0.0001 DZEO 2*H

NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "50000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.6.c -o $INIT -lm
./$INIT -ri 15. -ro 40 -rdzeo 24. -drdzeo 2. -sigma0 0.0001 -index 0.5 -alpha 0.01 -m0 1 -n 500

DRIFT=drift$NOW

gcc -O3 dustdrift1Dver14.c -o $DRIFT -lm
./$DRIFT -n 500 

wait

