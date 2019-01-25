#!/bin/bash

.  /local/env/envfastsimbac-1.0.1_bd3ad13d8f79.sh 

mkdir "FILE"$SLURM_ARRAY_TASK_ID
cp segsite "FILE"$SLURM_ARRAY_TASK_ID
cp "priorsAM-"$SLURM_ARRAY_TASK_ID  "FILE"$SLURM_ARRAY_TASK_ID
cd "FILE"$SLURM_ARRAY_TASK_ID

   t=`cat "priorsAM-"$SLURM_ARRAY_TASK_ID |cut -d" " -f1`
   r=`cat "priorsAM-"$SLURM_ARRAY_TASK_ID |cut -d" " -f2`
   delta=`cat "priorsAM-"$SLURM_ARRAY_TASK_ID |cut -d" " -f3`
   TS=`cat "priorsAM-"$SLURM_ARRAY_TASK_ID |cut -d" " -f4`
   N2=`cat "priorsAM-"$SLURM_ARRAY_TASK_ID |cut -d" " -f5`
   M12=`cat "priorsAM-"$SLURM_ARRAY_TASK_ID | cut -d" " -f7`
   M21=`cat "priorsAM-"$SLURM_ARRAY_TASK_ID | cut -d" " -f8`
   Tam=`cat "priorsAM-"$SLURM_ARRAY_TASK_ID | cut -d" " -f9`

echo   ./fastSimBac 30 1712736 -T -t $t -r $r $delta -I 2 13 17 0 -m 1 2 0 -m 2 1 0 -n  1 1  -n 2 $N2 -ema $Tam 2 0  $M12 $M21 0  -ej $TS 1 2     > $SLURM_ARRAY_TASK_ID".line"
   fastSimBac 30 1712736 -T -t $t -r $r $delta -I 2 13 17 0 -m 1 2 0 -m 2 1 0 -n  1 1  -n 2 $N2 -ema $Tam 2 0  $M12 $M21 0  -ej $TS 1 2     | msformatter > $SLURM_ARRAY_TASK_ID".temp2" 
   if test -s $SLURM_ARRAY_TASK_ID.temp2
   then echo  $SLURM_ARRAY_TASK_ID".temp2" >>log
   cat   $SLURM_ARRAY_TASK_ID".line" $SLURM_ARRAY_TASK_ID".temp2"  >> "priorsAM-"$SLURM_ARRAY_TASK_ID.ms
   else cat $SLURM_ARRAY_TASK_ID".line" segsite >>  "priorsAM-"$SLURM_ARRAY_TASK_ID.ms 
   fi
   cd ..
