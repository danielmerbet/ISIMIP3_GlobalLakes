#!/bin/bash

#SBATCH --job-name=ISIMIP3a-run
#SBATCH --output=%x-%j.out
#SBATCH --nodes=3
#SBATCH --ntasks=100
#SBATCH --time=05-00:00:00

echo EMPIEZA TODO `date`
for i in runs/*.sh;do 
  $i &
done
wait
echo TERMINA TODO `date`
