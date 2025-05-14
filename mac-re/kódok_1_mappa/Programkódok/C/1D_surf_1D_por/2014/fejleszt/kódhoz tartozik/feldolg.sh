#!/bin/bash

#sh abra.sh
#sh abra_20.sh

cat pormozgas.dat | grep 199800 > por_r.dat

gcc rand.c -o rand -lm
./rand

