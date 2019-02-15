#!/bin/bash
module load bioinfo/fastsimbac-bd3ad13
for i in `cat listems`
do
while read line 
  do 
  L=`echo $line |cut -d" " -f1`
  t=`echo $line |cut -d" " -f2`
  r=`echo $line |cut -d" " -f3`
  delta=`echo $line |cut -d" " -f4`
  TS=`echo $line |cut -d" " -f5`
  N2=`echo $line |cut -d" " -f6`
  Na=`echo $line |cut -d" " -f7`


    ./msnsam 30 $L -T -t $t  -I 2 13 17  -n 1 1 -n 2 $N2   -ej $TS 1 2  -eN $TS $Na > $i.ms
done <$i
done
