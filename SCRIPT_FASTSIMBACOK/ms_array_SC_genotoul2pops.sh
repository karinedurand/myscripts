#!/bin/bash

module load bioinfo/ms_2017
module load bioinfo/msums-ca90e3a

mkdir "FILE"$SLURM_ARRAY_TASK_ID
cp SC-$SLURM_ARRAY_TASK_ID "FILE"$SLURM_ARRAY_TASK_ID
cp spinput "FILE"$SLURM_ARRAY_TASK_ID
cd "FILE"$SLURM_ARRAY_TASK_ID
while read line 
do L=`echo $line |cut -d" " -f1`
t=`echo $line |cut -d" " -f2`
r=`echo $line |cut -d" " -f3`
#delta=`echo $line |cut -d" " -f4`
TS=`echo $line |cut -d" " -f4`
N1=`echo $line |cut -d" " -f5`
N2=`echo $line |cut -d" " -f6`
Na=`echo $line |cut -d" " -f7`
M12=`echo $line | cut -d" " -f8`
M21=`echo $line | cut -d" " -f9`
Tsc=`echo $line | cut -d" " -f10`
ms 30 1  -t $t -r $r  $L  -I 2 13 17   -n 1 $N1  -n 2 $N2   -m  1 2  $M12  -m  2 1  $M21 -eM $Tsc  0  -ej $TS 1 2  -eN $TS $Na >> SC-$SLURM_ARRAY_TASK_ID.ms
done< SC-$SLURM_ARRAY_TASK_ID 

echo SC-$SLURM_ARRAY_TASK_ID.ms > temp
cat spinput temp >spinput.txt
msums -S all -i spinput.txt -o ABCstat.txt
cd ..
