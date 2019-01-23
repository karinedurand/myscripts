#!/bin/bash
module load bioinfo/fastsimbac-bd3ad13
for i in `cat listems`
do
while read line 
  do   t=`echo $line |cut -d" " -f1`
  r=`echo $line |cut -d" " -f2`
  delta=`echo $line |cut -d" " -f3`
  TS=`echo $line |cut -d" " -f4`
  N2=`echo $line |cut -d" " -f5`
  Na=`echo $line |cut -d" " -f6`
  TS2=`echo $line |cut -d" " -f7`
  echo   ./fastSimBac 30 1712736 -T -t $t -r $r $delta -n 1 1 -n 2 $N2 -I 2 13 17  -ej $TS 1 2   > line
    ./fastSimBac 30 1712736 -T -t $t -r $r $delta -I 2 13 17  -n 1 1 -n 2 $N2   -ej $TS 1 2   | ./msformatter >temp2 
if test -s temp2
then 
echo "temp2" >>log
cat  line temp2  >> $i.ms
else
cat line segsite >>  $i.ms
fi
done< $i
done