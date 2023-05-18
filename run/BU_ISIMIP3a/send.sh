#!/bin/bash
#SBATCH --job-name=ISIMIP3a-run
#SBATCH --output=%x-%j.out
#SBATCH --ntasks=100
#SBATCH --time=05-00:00:00

#module load parallel
ml parallel/20210722-GCCcore-11.2.0

echo EMPIEZA TODO `date`
echo 20CRv3-W5E5 obsclim
parallel -j $SLURM_NTASKS srun -N 1 -n 1 -c 1 --exact bash ::: runs/*.sh

echo TERMINA TODO `date`
