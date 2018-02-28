#!/bin/bash
#SBATCH --ntasks=136
#SBATCH --ntasks-per-node=68
#SBATCH --time=24:00:00
#SBATCH --account=Pra15_3924
#SBATCH -p knl_usr_prod
#SBATCH -e para-RT.e
#SBATCH -o para-RT.o

cd $SLURM_SUBMIT_DIR

BASE_DIR=$CINECA_SCRATCH"/reditools/"
INPUT_DIR="/marconi_scratch/userexternal/epicardi/PRJNA231202/SRR1047874/"
OUTPUT_DIR=$BASE_DIR"/output/"

SAMPLE_ID="SRR1047874"
SOURCE_BAM_FILE=$INPUT_DIR$SAMPLE_ID".bam"

REFERENCE=$BASE_DIR"hg19.fa"
OMOPOLYMER_FILE=$BASE_DIR"omopolymeric_positions.txt"
SIZE_FILE=$BASE_DIR"hg19.chrom.sizes"

COVERAGE_DIR=$BASE_DIR"/cov/"$SAMPLE_ID"/"
COVERAGE_FILE=$COVERAGE_DIR$SAMPLE_ID".cov"
TEMP_DIR=$BASE_DIR"/temp/"$SAMPLE_ID"/"

OUTPUT=$OUTPUT_DIR/$SAMPLE_ID/table.gz
NUM_CORES=68

echo "Launching REDItool on $SAMPLE_ID (output_dir=$OUTPUT_DIR)";
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

if [ ! -f $COVERAGE_FILE ]
then
        t1=$(date +%s)
        echo "[STATS] [COVERAGE] START="$t1
        ./extract_coverage.sh $SOURCE_BAM_FILE $COVERAGE_DIR $SIZE_FILE
        t2=$(date +%s)
        elapsed_time=$(($t2-$t1))
        elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
        echo "[STATS] [COVERAGE] START="$t1" END="$t2" ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human
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
time mpirun src/cineca/parallel_reditools.py -f $SOURCE_BAM_FILE -r $REFERENCE -m $OMOPOLYMER_FILE -G $COVERAGE_FILE -D $COVERAGE_DIR -t $TEMP_DIR -Z $SIZE_FILE $options 2>&1 | tee $SAMPLE_ID.log
t2=$(date +%s)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [PARALLEL] START="$t1" END="$t2" ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human

t1=$(date +%s)
time ./merge.sh $TEMP_DIR $OUTPUT $NUM_CORES
t2=$(date +%s)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [MERGE] START="$t1" END="$t2" ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human

echo "END:"`date`
