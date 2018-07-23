module load autoload profile/global
module load ig_homo_sapiens/hg19
REFERENCE=$IG_HG19_GENOME"/genome.fa"
SIZE_FILE=$REFERENCE".fai"

SAMPLE_FILE=$BASE_DIR"samples-10.txt"

N=$(cat $SAMPLE_FILE | wc -l | cut -d' ' -f 1)

BASE_DIR="/marconi_scratch/userexternal/tflati00/reditools_paper/"
COVERAGE_DIR=$BASE_DIR"/cov-multisample/"

for SOURCE_BAM_FILE in $(cat $SAMPLE_FILE | head -n 1)
do
	if [ ! -s $SOURCE_BAM_FILE ]
	then
		echo "File $SOURCE_BAM_FILE does not exists. Skipping."
		continue
	fi
	
	SAMPLE_ID=$(basename $SOURCE_BAM_FILE | sed 's/\.bam//g')
	COV=$COVERAGE_DIR$SAMPLE_ID"/"
	COV_FILE=$COV$SAMPLE_ID".cov"

	date

	if [ ! -f $COV_FILE ]
	then
			echo "[STATS] [COVERAGE] [$SAMPLE_ID]"
			sbatch --export=ALL -o $COV/output.txt -e $COV/error.txt reditools2.0/extract_coverage_slurm.sh &
	fi
done