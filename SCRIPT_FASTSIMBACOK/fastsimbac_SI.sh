#!/bin/bash
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
  echo   ./fastSimBac 30 $L -T -t $t -r $r $delta -I 2 13 17  -ej $TS 1 2  -eN $TS $Na > line
    ./fastSimBac 30 $L -T -t $t -r $r $delta -I 2 13 17  -ej $TS 1 2  -eN $TS $Na | ./msformatter >temp2 
if test -s temp2
then 
echo "temp2" >>log
cat  line temp2  >> $i.ms
else
cat line segsite >>  $i.ms
fi
done< $i
done