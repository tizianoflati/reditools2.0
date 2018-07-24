#!/bin/bash
#SBATCH --ntasks=136
#SBATCH --ntasks-per-node=68
#SBATCH --time=4:00:00
#SBATCH --account=Pra15_3924
#SBATCH -p knl_usr_prod
#SBATCH -e para-RT-MT.e
#SBATCH -o para-RT-MT.o

cd $SLURM_SUBMIT_DIR

BASE_DIR="/marconi_scratch/userexternal/tflati00/reditools_paper/"

OUTPUT_DIR=$BASE_DIR"/output-multisample-2/"
TEMP_DIR=$BASE_DIR"/tmp-multisample-2/"
COVERAGE_DIR=$BASE_DIR"/cov-multisample-2/"

#DATA_DIR="/home/flati/data/reditools/input/"
DATA_DIR="$CINECA_SCRATCH/public/"

module load autoload profile/global
module load ig_homo_sapiens/hg19
REFERENCE=$IG_HG19_GENOME"/genome.fa"
#REFERENCE=$DATA_DIR"hg19m.fa"

OMOPOLYMER_FILE=$DATA_DIR"omopolymeric_positions.txt"
SIZE_FILE=$REFERENCE".fai"

SAMPLE_FILE=$BASE_DIR"samples.txt"

# NUM_CORES=68

if [ ! -s $SAMPLE_FILE ]
then
	echo "File $SAMPLE_FILE does not exist. Please, provide an existing file."
	exit
fi

# Environment setup
module load python/2.7.12
source ENV/bin/activate
module load autoload openmpi/1-10.3--gnu--6.1.0
# module load autoload samtools
module load autoload htslib

# for SOURCE_BAM_FILE in $(cat $SAMPLE_FILE)
# do
# 	if [ ! -s $SOURCE_BAM_FILE ]
# 	then
# 		echo "File $SOURCE_BAM_FILE does not exists. Skipping."
# 		continue
# 	fi
# 	
# 	SAMPLE_ID=$(basename $SOURCE_BAM_FILE | sed 's/\.bam//g')
# 	COV=$COVERAGE_DIR$SAMPLE_ID"/"
# 	COV_FILE=$COV$SAMPLE_ID".cov"
# 
# 	date
# 
# 	if [ ! -f $COV_FILE ]
# 	then
# 			echo "Launching REDItool COVERAGE on $SAMPLE_ID (output_dir=$COV)";
# 			
# 			t1=$(date +%s)
# 			t1_human=$(date)
# 			echo "[STATS] [COVERAGE] [$SAMPLE_ID] START="$t1_human" ["$t1"]"
# 			time ./extract_coverage.sh $SOURCE_BAM_FILE $COV $SIZE_FILE &
# 			t2=$(date +%s)
# 			t2_human=$(date)
# 			elapsed_time=$(($t2-$t1))
# 			elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
# 			echo "[STATS] [COVERAGE] [$SAMPLE_ID] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human
# 	fi
# done
# wait

# strand=0
# options=""
# if [ $strand != 0 ]
# then
#         options="-C -T 2 -s $strand"
# fi
options=""

COV_FILE=$COV$SAMPLE_ID".cov"
TEMP=$TEMP_DIR$SAMPLE_ID"/"
OUTPUT=$OUTPUT_DIR/$SAMPLE_ID/table.gz

# Program launch
echo "START:"`date`
t1=$(date +%s)
t1_human=$(date)
# time mpirun src/cineca/reditools2_multisample.py -F $SAMPLE_FILE -r $REFERENCE -m $OMOPOLYMER_FILE -D $COVERAGE_DIR -t $TEMP_DIR -Z $SIZE_FILE $options 2>&1 | tee MULTI_SAMPLES.log
time mpirun src/cineca/reditools2_multisample.py -F $SAMPLE_FILE -r $REFERENCE -D $COVERAGE_DIR -t $TEMP_DIR -Z $SIZE_FILE $options 2>&1 | tee MULTI_SAMPLES.log
t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [PARALLEL] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human

# export PATH=$HTSLIB_HOME/bin/:$PATH
# for SOURCE_BAM_FILE in $(cat $SAMPLE_FILE)
# do
# 	t1=$(date +%s)
# 	t1_human=$(date)
# 	
# 	SAMPLE_ID=$(basename $SOURCE_BAM_FILE | sed 's/\.bam//g')
# 
# 	COV=$COVERAGE_DIR$SAMPLE_ID"/"	
# 	COV_FILE=$COV$SAMPLE_ID".cov"
# 	TEMP=$TEMP_DIR$SAMPLE_ID"/"
# 	OUTPUT=$OUTPUT_DIR/$SAMPLE_ID/table.gz
# 
# 	time ./merge.sh $TEMP $OUTPUT $NUM_CORES &
# 	t2=$(date +%s)
# 	t2_human=$(date)
# 	elapsed_time=$(($t2-$t1))
# 	elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
# 	echo "[STATS] [MERGE] [$SAMPLE_ID] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human
# 
# 	echo "[$SAMPLE_ID] END:"`date`
# 	echo "OK" > $TEMP/status.txt
# done
# wait
