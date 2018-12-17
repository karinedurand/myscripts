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
  echo "l" $L
  echo "t" $t
  echo "r"  $r
  echo "delta" $delta
  echo "N1"  $N1
  echo "n2" $N2
  echo "TS" $TS
  echo   ./fastSimBac 30 $L -T -t $t -r $r $delta -I 2 13 17  -ej $TS 1 2  > line
    ./fastSimBac 30 $L -T -t $t -r $r $delta -I 2 13 17  -ej $TS 1 2  | ./msformatter >temp2 
if test -s temp2
then 
echo "temp2" >>log
cat  line temp2  >> $i.ms
else
cat line segsite >>  $i.ms
fi
done< $i
done