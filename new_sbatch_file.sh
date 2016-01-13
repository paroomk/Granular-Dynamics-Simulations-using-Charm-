#!/bin/sh
#SBATCH --partition=slurm_me759
#SBATCH --time=0-00:19:59
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --output=./leanmd.o%j
./charmrun +p16 ./leanmd  16 16 32 32 2000 output5
