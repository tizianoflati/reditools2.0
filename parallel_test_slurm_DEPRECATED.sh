#!/bin/bash

#SBATCH --job-name=REDItools2Job
#SBATCH -N 1
#SBATCH -n 36
#SBATCH -p gll_usr_prod
#SBATCH --mem=115GB
#SBATCH --time 05:00:00
#SBATCH --account ELIX4_manniron
#SBATCH --error REDItools2Job.err
#SBATCH --output REDItools2Job.out


#########################################################
########    Parameters setting
#########################################################

##SAMPLE_ID is the basename of the sample of interest
SAMPLE_ID="SRR2135332"

##bam file to be analysed
SOURCE_BAM_FILE="test/SRR2135332.bam"

##reference chromosome or genome
REFERENCE="test/chr21.fa"
REFERENCE_DNA=$(basename "$REFERENCE")

##fasta index file created by samtools
SIZE_FILE="test/chr21.fa.fai"

##number of utilized cores
NUM_CORES=2

##setting output file
OUTPUT_FILE="test_results/output/parallel_table.txt.gz"
TEMP_DIR="test_results/temp/"

##setting the coverage file
COVERAGE_FILE="test_results/coverage/SRR2135332.cov"

##setting coverage directory
COVERAGE_DIR="test_results/coverage/"

##setting output directory
OUTPUT_DIR=$(basename "$OUTPUT_FILE")

#########################################################
########    Modules loading
#########################################################

module load profile/bioinf
module load python/2.7.12
module load autoload samtools/1.9
module load autoload profile/global
module load autoload openmpi/3.1.4--gnu--7.3.0
module load autoload samtools

echo "Launching REDItool on $SAMPLE_ID (output_file=$OUTPUT_FILE)";

#########################################################
########    Coverage
#########################################################

## If the coverage file doesn’t exist, then the script calculate It.

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

#########################################################
########    Parallel Computation
#########################################################

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

time mpirun src/cineca/parallel_reditools.py -g $REFERENCE_DNA -f $SOURCE_BAM_FILE -r $REFERENCE -G $COVERAGE_FILE -D $COVERAGE_DIR -t $TEMP_DIR -Z $SIZE_FILE $options 2>&1 | tee $SAMPLE_ID.log
t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [PARALLEL] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human

#########################################################
########    Merging
#########################################################

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
