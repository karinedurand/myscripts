#!/bin/bash

.  /local/env/envfastsimbac-1.0.1_bd3ad13d8f79.sh 

mkdir "FILE"$SLURM_ARRAY_TASK_ID
cp segsite "FILE"$SLURM_ARRAY_TASK_ID
cp "PriorsSC-"$SLURM_ARRAY_TASK_ID  "FILE"$SLURM_ARRAY_TASK_ID
cd "FILE"$SLURM_ARRAY_TASK_ID

   t=`cat "PriorsSC-"$SLURM_ARRAY_TASK_ID |cut -d" " -f1`
   r=`cat "PriorsSC-"$SLURM_ARRAY_TASK_ID |cut -d" " -f2`
   delta=`cat "PriorsSC-"$SLURM_ARRAY_TASK_ID |cut -d" " -f3`
   TS=`cat "PriorsSC-"$SLURM_ARRAY_TASK_ID |cut -d" " -f4`
   N2=`cat "PriorsSC-"$SLURM_ARRAY_TASK_ID |cut -d" " -f5`
   M12=`cat "PriorsSC-"$SLURM_ARRAY_TASK_ID | cut -d" " -f7`
   M21=`cat "PriorsSC-"$SLURM_ARRAY_TASK_ID | cut -d" " -f8`
   Tsc=`cat "PriorsSC-"$SLURM_ARRAY_TASK_ID  | cut -d" " -f9`
   echo   ./fastSimBac 30 1712736  -T -t $t -r $r $delta -I 2 13 17   -n 1 1 -n 2 $N2   -m  1 2  $M12  -m  2 1  $M21 -eM $Tsc  0  -ej $TS 1 2    > $SLURM_ARRAY_TASK_ID".line"
   fastSimBac 30 1712736 -T -t $t -r $r $delta -I 2 13 17   -n 1 1 -n 2 $N2   -m  1 2  $M12  -m  2 1  $M21 -eM $Tsc  0  -ej $TS 1 2    | msformatter > $SLURM_ARRAY_TASK_ID".temp2" 
   if test -s $SLURM_ARRAY_TASK_ID.temp2
   then echo  $SLURM_ARRAY_TASK_ID".temp2" >>log
   cat   $SLURM_ARRAY_TASK_ID".line" $SLURM_ARRAY_TASK_ID".temp2"  >> "PriorsSC-"$SLURM_ARRAY_TASK_ID.ms
   else cat $SLURM_ARRAY_TASK_ID".line" segsite >>  "PriorsSC-"$SLURM_ARRAY_TASK_ID.ms 
   fi
   cd ..
