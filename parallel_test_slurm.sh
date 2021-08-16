#!/bin/bash
#SBATCH --job-name=REDItools2Job
#SBATCH -N 3
#SBATCH -n 12
#SBATCH -p m100_usr_prod
#SBATCH --time 02:00:00
#SBATCH --account cin_staff
#SBATCH --error REDItools2Job.err
#SBATCH --output REDItools2Job.out

SAMPLE_ID="SRR2135332"
SOURCE_BAM_FILE="test/SRR2135332.bam"
REFERENCE="test/chr21.fa"
REFERENCE_DNA=$(basename "$REFERENCE")
SIZE_FILE="test/chr21.fa.fai"
NUM_CORES=12
OUTPUT_FILE="test_results/output/parallel_table.txt.gz"
TEMP_DIR="test_results/temp/"
COVERAGE_FILE="test_results/coverage/SRR2135332.cov"
COVERAGE_DIR="test_results/coverage/"
OUTPUT_DIR=$(basename "$OUTPUT_FILE")


module load spack
module load python/2.7.16--gcc--8.4.0-bgv
module load autoload py-mpi4py/3.0.3--gcc--8.4.0-spectrmpi-ac2
module load py-virtualenv/16.7.6--gcc--8.4.0-4ut
module load profile/global
module load samtools/1.12

source ENV/bin/activate

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

time mpirun -np $NUM_CORES src/cineca/parallel_reditools.py -f $SOURCE_BAM_FILE -o $OUTPUT_FILE -r $REFERENCE -t $TEMP_DIR -Z $SIZE_FILE -G $COVERAGE_FILE -D $COVERAGE_DIR $options 2>&1 | tee $SAMPLE_ID.log

t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [PARALLEL] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human

t1=$(date +%s)
t1_human=$(date)
export PATH=$HTSLIB_HOME/bin/:$PATH
time ./merge.sh $TEMP_DIR $OUTPUT_FILE $NUM_CORES
t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [MERGE] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human

echo "END:"`date`
echo "OK" > $TEMP_DIR/status.txt

deactivate
