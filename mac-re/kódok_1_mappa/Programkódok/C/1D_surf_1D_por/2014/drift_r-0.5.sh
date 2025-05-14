#!/bin/bash

TIME=time.dat

# ALPHA 0.0001 DZEO 2*H

NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "8000. 150. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 7. -sigma0 0.000135013 -index -0.5 -rdzeo 5. -drdzeo 1.0 -flind 0.0 -alpha 0.01 -m0 2 -n 500

DRIFT=drift$NOW

gcc -O3 dustdrift1Dver12.c -o $DRIFT -lm
#nohup ./$DRIFT &

wait

