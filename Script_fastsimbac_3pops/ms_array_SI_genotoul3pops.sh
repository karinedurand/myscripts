#!/bin/bash

module load bioinfo/ms_2017

while read line 
  do  L=`echo $line |cut -d" " -f1`
  t=`echo $line |cut -d" " -f2`
  r=`echo $line |cut -d" " -f3`
  TS1=`echo $line |cut -d" " -f4`
  TS2=`echo $line |cut -d" " -f5`
  N1=`echo $line |cut -d" " -f6`
  N2=`echo $line |cut -d" " -f7`
  N3=`echo $line |cut -d" " -f8`
  Na=`echo $line |cut -d" " -f9`
  
  
    ms 36 1  -t $t -r $r $L  -I 3 17 13 6 0  -n 1 $N1 -n 2 $N2 -n 3 $N3  -ej $TS2 2 3 -ej $TS1 3 1 -eN $TS2 $Na  >>$SLURM_ARRAY_TASK_ID.ms
    
done< SI_A-$SLURM_ARRAY_TASK_ID 
