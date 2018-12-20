#!/bin/bash
module load bioinfo/fastsimbac-bd3ad13
####SC   : 
#L locus lengh
#t theta
#r rho
# λ length of recombining intervals
#N1 size pop1
#N2 size pop2

#Ts time scale of split SC
for i in `cat listems`
do while read line 
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
Tsc=`echo $line | cut -d" " -f11`
TS2=`echo $line | cut -d" " -f12`
echo ./fastSimBac 30 $L -T -t $t -r $r $delta -I 2 13 17   -n 1 $N1 -n 2 $N2   -m  1 2  $M12  -m  2 1  $M21 -eM $Tsc  0  -ej $TS 1 2 -eN $TS2 $Na  > line
fastSimBac 30 $L  -t $t -r $r $delta -I 2 13 17   -n 1 $N1 -n 2 $N2   -m  1 2  $M12  -m  2 1  $M21 -eM $Tsc  0  -ej $TS 1 2  -eN $TS2 $Na | msformatter >temp2 
if test -s temp2
then 
echo "temp2" >>log
cat  line temp2  >> $i.ms
else
cat line segsite >>  $i.ms
fi
done< $i
done
