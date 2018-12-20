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
  TS2=`echo $line |cut -d" " -f9`
  echo   ./fastSimBac 30 $L -T -t $t -r $r $delta -I 2 13 17  -ej $TS 1 2  -eN $TS2 $Na > line
    ./fastSimBac 30 $L -T -t $t -r $r $delta -I 2 13 17  -ej $TS 1 2  -eN $TS2 $Na | ./msformatter >temp2 
if test -s temp2
then 
echo "temp2" >>log
cat  line temp2  >> $i.ms
else
cat line segsite >>  $i.ms
fi
done< $i
done