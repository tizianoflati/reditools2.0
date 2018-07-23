#!/bin/bash
#SBATCH --ntasks=68
#SBATCH --ntasks-per-node=68
#SBATCH --time=02:00:00
##SBATCH --account=Pra15_3924
#SBATCH --account=cin_staff
#SBATCH -p knl_usr_prod

# SAMPLE_ID
# OUTPUT_DIR
# SOURCE_BAM_FILE
# COV
# SIZE_FILE

cd $SLURM_SUBMIT_DIR

echo "Launching REDItool COVERAGE on $SAMPLE_ID (output_dir=$OUTPUT_DIR)";

t1=$(date +%s)
t1_human=$(date)
echo "[STATS] [COVERAGE] [$SAMPLE_ID] START="$t1_human" ["$t1"]"
time ./extract_coverage.sh $SOURCE_BAM_FILE $COV $SIZE_FILE &
t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [COVERAGE] [$SAMPLE_ID] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human