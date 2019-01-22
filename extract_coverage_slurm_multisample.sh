module load autoload profile/global
module load ig_homo_sapiens/hg19
REFERENCE=$IG_HG19_GENOME"/genome.fa"
export SIZE_FILE=$REFERENCE".fai"

BASE_DIR="/marconi_scratch/userexternal/tflati00/reditools_paper/"
SAMPLE_FILE=$BASE_DIR"samples-10.txt"
N=$(cat $SAMPLE_FILE | wc -l | cut -d' ' -f 1)

COVERAGE_DIR=$BASE_DIR"/cov-multisample/"

for SOURCE_BAM_FILE in $(cat $SAMPLE_FILE | head -n 1)
do
	if [ ! -s $SOURCE_BAM_FILE ]
	then
		echo "File $SOURCE_BAM_FILE does not exists. Skipping."
		continue
	fi
	
	export SOURCE_BAM_FILE
	export SAMPLE_ID=$(basename $SOURCE_BAM_FILE | sed 's/\.bam//g')
	export COV=$COVERAGE_DIR$SAMPLE_ID"/"
	export COV_FILE=$COV$SAMPLE_ID".cov"

	mkdir -p $COV

	if [ ! -f $COV_FILE ]
	then
			echo "[STATS] [COVERAGE] [$SAMPLE_ID]"
			sbatch --export=ALL -J cov-$SAMPLE_ID -o $COV/output.txt -e $COV/error.txt ./extract_coverage_slurm.sh &
	fi
done