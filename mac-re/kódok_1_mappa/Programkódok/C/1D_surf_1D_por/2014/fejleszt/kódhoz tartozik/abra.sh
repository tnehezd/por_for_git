#!/bin/bash


POR="pormozgas.dat"
TEMP="temp1.dat"
TEMP2="temp2.dat"
GPLOT="plot.por"
OUTPUT="kesz.dat"


touch $TEMP
touch $TEMP2

cat $POR | awk '{print $2,$3,$1 }' | grep -v '250$' | grep -v '750$' > $TEMP
sort -g -k 1 -k 3 $TEMP > $TEMP2
rm $TEMP

gcc break.c
./a.out $OUTPUT

rm $TEMP2

rm $GPLOT
touch $GPLOT

echo "set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )" >> $GPLOT
echo "set term png size 1200, 800" >> $GPLOT
echo "set output \"abra.png\"" >> $GPLOT
echo "set xlabel \"Time (10^3 yr)\"" >> $GPLOT
echo "set ylabel \"R (au)\"" >> $GPLOT
echo "plot \"$OUTPUT\" every :::0::1 u (\$3/1000):2:1 w l palette title \"\"" >> $GPLOT

echo "load \"$GPLOT\"" | gnuplot
