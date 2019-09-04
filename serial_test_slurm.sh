#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --account=cin_staff
#SBATCH -p knl_usr_prod
#SBATCH -e serial-RT.e
#SBATCH -o serial-RT.o

# Serial test (SLURM)#
module load python/2.7.12

source ENV/bin/activate

python src/cineca/reditools.py -f test/SRR2135332.bam -g chr21 -r test/chr21.fa -o serial_table_slurm.txt

deactivate
