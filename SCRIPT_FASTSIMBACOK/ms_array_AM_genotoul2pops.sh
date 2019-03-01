#!/bin/bash

module load bioinfo/ms_2017

mkdir "FILE"$SLURM_ARRAY_TASK_ID
cp AM-$SLURM_ARRAY_TASK_ID "FILE"$SLURM_ARRAY_TASK_ID
cd "FILE"$SLURM_ARRAY_TASK_ID
while read line 
  do 
  L=`echo $line |cut -d" " -f1`
  t=`echo $line |cut -d" " -f2`
  r=`echo $line |cut -d" " -f3`
  TS=`echo $line |cut -d" " -f4`
  N1=`echo $line |cut -d" " -f5`
  N2=`echo $line |cut -d" " -f6`
  Na=`echo $line |cut -d" " -f7`
  M12=`echo $line |cut -d" " -f8`
  M21=`echo $line |cut -d" " -f9`
  Tam=`echo $line |cut -d" " -f10`

    ms 30 1  -t $t -r $r $L  -I 2 13 17  -n 1 $N1 -n 2 $N2   -m 1 2 0 -m 2 1 0  -ema $Tam 2 0  $M12 $M21 0  -ej $TS 1 2  -eN $TS $Na   >> AM-$SLURM_ARRAY_TASK_ID.ms
   done< AM-$SLURM_ARRAY_TASK_ID 
cd ..
