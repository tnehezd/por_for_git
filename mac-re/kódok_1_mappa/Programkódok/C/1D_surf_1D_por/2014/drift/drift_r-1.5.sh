#!/bin/bash


TIME=time.dat

# ALPHA 0.0001 DZEO 2*H

NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 2.0 -flind 0.25 -alpha 0.0001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 2.0 -flind 0.125 -alpha 0.0001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 2.0 -flind 0.0 -alpha 0.0001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


# ALPHA 0.0001 DZEO 1.5*H

NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 1.5 -flind 0.25 -alpha 0.0001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 1.5 -flind 0.125 -alpha 0.0001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 1.5 -flind 0.0 -alpha 0.0001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait



# ALPHA 0.001 DZEO 2*H

NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 2.0 -flind 0.25 -alpha 0.001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 2.0 -flind 0.125 -alpha 0.001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 2.0 -flind 0.0 -alpha 0.001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


# ALPHA 0.001 DZEO 1.5*H

NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 1.5 -flind 0.25 -alpha 0.001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 1.5 -flind 0.125 -alpha 0.001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 1.5 -flind 0.0 -alpha 0.001 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait



# ALPHA 0.01 DZEO 2*H

NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 2.0 -flind 0.25 -alpha 0.01 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 2.0 -flind 0.125 -alpha 0.01 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 2.0 -flind 0.0 -alpha 0.01 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


# ALPHA 0.01 DZEO 1.5*H

NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 1.5 -flind 0.25 -alpha 0.01 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 1.5 -flind 0.125 -alpha 0.01 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait


NOW=$(date +"%m%d%H%M%S")
touch $TIME
echo "300000. 250. $NOW" > $TIME

INIT=init$NOW

gcc init_dust.ver.5.c -o $INIT -lm
./$INIT -ri 2 -ro 150 -sigma0 0.000135013 -index -1.5 -rdzeo 35 -drdzeo 1.5 -flind 0.0 -alpha 0.01 -m0 2 -n 1000

DRIFT=drift$NOW

echo $DRIFT


gcc -O3 dustdrift1Dver11.c -o $DRIFT -lm
nohup ./$DRIFT &

wait

