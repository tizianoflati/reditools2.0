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
        ./extract_coverage.sh $SOURCE_BAM_FILE $COVERAGE_DIR $SIZE_FILE
fi

strand=0
options=""
if [ $strand != 0 ]
then
        options="-C -T 2 -s $strand"
fi

# Program launch
echo "START:"`date`

time mpirun src/cineca/parallel_reditools.py -f $SOURCE_BAM_FILE -r $REFERENCE -m $OMOPOLYMER_FILE -o $OUTPUT_DIR/$SAMPLE_ID/table.gz -G $COVERAGE_FILE -D $COVERAGE_DIR -t $TEMP_DIR -Z $SIZE_FILE $options 2>&1 | tee $SAMPLE_ID.log

echo "END:"`date`
