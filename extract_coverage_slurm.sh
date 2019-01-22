#!/bin/bash
#SBATCH --ntasks=25
#SBATCH --ntasks-per-node=25
#SBATCH --time=02:00:00
##SBATCH --account=Pra15_3924
#SBATCH --account=cin_staff
#SBATCH -p knl_usr_prod

# SAMPLE_ID
# SOURCE_BAM_FILE
# COV
# SIZE_FILE

echo "Launching REDItool COVERAGE on $SAMPLE_ID";

module load autoload profile/global
module load autoload samtools

t1=$(date +%s)
t1_human=$(date)
echo "[STATS] [COVERAGE] [$SAMPLE_ID] START="$t1_human" ["$t1"]"
time ./extract_coverage_dynamic.sh $SOURCE_BAM_FILE $COV $SIZE_FILE
t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [COVERAGE] [$SAMPLE_ID] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human
