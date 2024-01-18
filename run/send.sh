#!/bin/bash
#SBATCH --job-name=ISIMIP3bR
#SBATCH --output=%x-%j.out
#SBATCH --ntasks=100
#SBATCH --time=05-00:00:00


###SBATCH --partition=broadwell

#module load parallel
ml parallel/20210722-GCCcore-11.2.0

echo EMPIEZA TODO `date`
parallel -j $SLURM_NTASKS srun -N 1 -n 1 -c 1 --exact bash ::: runs/*.sh

echo TERMINA TODO `date`
