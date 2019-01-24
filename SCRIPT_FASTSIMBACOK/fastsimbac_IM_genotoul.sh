---#!/bin/bash

module load bioinfo/fastsimbac-bd3ad13


mkdir "FILE"$SLURM_ARRAY_TASK_ID
cp segsite "FILE"$SLURM_ARRAY_TASK_ID
mv "priorsIM-"$SLURM_ARRAY_TASK_ID  "FILE"$SLURM_ARRAY_TASK_ID
cd "FILE"$SLURM_ARRAY_TASK_ID

  t=`cat "priorsIM-"$SLURM_ARRAY_TASK_ID  |cut -d" " -f1`
  r=`cat "priorsIM-"$SLURM_ARRAY_TASK_ID |cut -d" " -f2`
  delta=`cat "priorsIM-"$SLURM_ARRAY_TASK_ID |cut -d" " -f3`
  TS=`cat "priorsIM-"$SLURM_ARRAY_TASK_ID |cut -d" " -f4`
  N2=`cat "priorsIM-"$SLURM_ARRAY_TASK_ID |cut -d" " -f5`
  M12=`cat "priorsIM-"$SLURM_ARRAY_TASK_ID | cut -d" " -f7`
  M21=`cat "priorsIM-"$SLURM_ARRAY_TASK_ID | cut -d" " -f8`

  echo   fastSimBac 30 1712736 -T -t $t -r $r $delta -n 1 1 -n 2 $N2 -I 2 13 17 -m  1 2  $M12  -m  2 1  $M21  -ej $TS 1 2   >  $SLURM_ARRAY_TASK_ID".line"
    fastSimBac 30 1712736 -T -t $t -r $r $delta -I 2 13 17  -n 1 1 -n 2 $N2  -m  1 2  $M12  -m  2 1  $M21  -ej $TS 1 2   | msformatter > $SLURM_ARRAY_TASK_ID".temp2"
if test -s $SLURM_ARRAY_TASK_ID.temp2
 then echo  $SLURM_ARRAY_TASK_ID".temp2" >>log
 cat   $SLURM_ARRAY_TASK_ID".line" $SLURM_ARRAY_TASK_ID".temp2"  >> $SLURM_ARRAY_TASK_ID.ms
 else cat $SLURM_ARRAY_TASK_ID".line" segsite >>  $SLURM_ARRAY_TASK_ID.ms 
 fi