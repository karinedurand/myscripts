#!/bin/bash
module load bioinfo/fastsimbac-bd3ad13
for i in `cat listems`
do
while read line 
  do L=`echo $line |cut -d" " -f1`
  t=`echo $line |cut -d" " -f2`
  r=`echo $line |cut -d" " -f3`
  delta=`echo $line |cut -d" " -f4`
  TS=`echo $line |cut -d" " -f5`
  N1=`echo $line |cut -d" " -f6`
  N2=`echo $line |cut -d" " -f7`
  Na=`echo $line |cut -d" " -f8`
  M12=`echo $line | cut -d" " -f9`
  M21=`echo $line | cut -d" " -f10`
  Tam=`echo $line | cut -d" " -f11`
  TS2=`echo $line | cut -d" " -f12`

./msnsam 30 $L  -t $t -r $r $delta -I 2 13 17 0 -m 1 2 0 -m 2 1 0 -n  1 $N1  -n 2 $N2 -ema $Tam 2 0  $M12 $M21 0  -ej $TS 1 2  -eN $TS2 $Na 
done< $i
done
