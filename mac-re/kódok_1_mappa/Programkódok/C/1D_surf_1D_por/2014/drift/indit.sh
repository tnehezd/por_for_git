#!/bin/bash

DIR1=R-0.5
DIR2=R-1.5
DIR3=R-1.0

mkdir $DIR1
mkdir $DIR2
mkdir $DIR3

cp *.c ./$DIR1
cp *0.5*.sh ./$DIR1
cp *.c ./$DIR2
cp *1.0*.sh ./$DIR2
cp *.c ./$DIR3
cp *1.5*.sh ./$DIR3

cd $DIR1
nohup sh *.sh &

cd .. 

cd $DIR2
nohup sh *.sh &

cd .. 

cd $DIR3
nohup sh *.sh &

FIN=README.dat

touch $FIN
echo "A FUTASOK ELINDULTAK!" > $FIN
