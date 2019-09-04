#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:20:00
#SBATCH --account=cin_staff
#SBATCH -p knl_usr_prod
#SBATCH -e para-RT.e
#SBATCH -o para-RT.o

##########################
### PARAMETER SETTING  ###
##########################
SAMPLE_ID="SRR2135332"
SOURCE_BAM_FILE="test/SRR2135332.bam"
REFERENCE="test/chr21.fa"
SIZE_FILE="test/chr21.fa.fai"

NUM_CORES=2
OUTPUT_FILE="test_results/output/parallel_table.txt.gz"
TEMP_DIR="test_results/temp/"
COVERAGE_FILE="test_results/coverage/SRR2135332.chr21.cov"
COVERAGE_DIR="test_results/coverage/"
OUTPUT_DIR=$(basename "$OUTPUT_FILE")

module load profile/global

echo "Launching REDItool on $SAMPLE_ID (output_file=$OUTPUT_FILE)";
date

if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir "$OUTPUT_DIR"
fi


# Environment setup
module load python/2.7.12
source ENV/bin/activate
module load autoload profile/global
module load autoload openmpi/1-10.3--gnu--6.1.0
module load autoload samtools
module load autoload htslib

##########################
### COVERAGE ANALYSIS  ###
##########################
if [ ! -f $COVERAGE_FILE ]
then
        t1=$(date +%s)
        t1_human=$(date)
        echo "[STATS] [COVERAGE] START="$t1_human" ["$t1"]"
        ./extract_coverage_dynamic.sh $SOURCE_BAM_FILE $COVERAGE_DIR $SIZE_FILE
        t2=$(date +%s)
        t2_human=$(date)
        elapsed_time=$(($t2-$t1))
        elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
        echo "[STATS] [COVERAGE] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human
fi


############################
### PARALLEL COMPUTATION ###
############################
strand=0
options=""
if [ $strand != 0 ]
then
        options="-C -T 2 -s $strand"
fi

# Program launch
echo "START:"`date`
t1=$(date +%s)
t1_human=$(date)
time mpirun src/cineca/parallel_reditools.py -g chr21 -f $SOURCE_BAM_FILE -r $REFERENCE -G $COVERAGE_FILE -D $COVERAGE_DIR -t $TEMP_DIR -Z $SIZE_FILE $options 2>&1 | tee $SAMPLE_ID.log
t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [PARALLEL] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human


######################
####### MERGE  #######
######################
t1=$(date +%s)
t1_human=$(date)
export PATH=$HTSLIB_HOME/bin/:$PATH
time ./merge.sh $TEMP_DIR $OUTPUT $NUM_CORES
t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [MERGE] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human

echo "END:"`date`
echo "OK" > $TEMP_DIR/status.txt

deactivate
