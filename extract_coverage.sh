COVERAGE_DIR=$2
FILENAME=`basename $1`
COVERAGE_FILE=${FILENAME%.*}

COVERAGE_DIR=$COVERAGE_DIR$COVERAGE_FILE"/"

if [ ! -d $COVERAGE_DIR ]
then
	mkdir -p $COVERAGE_DIR
fi

samtools depth $1 | grep -vP "\t0$" | tee $COVERAGE_DIR$COVERAGE_FILE".cov" | awk '{print($0) > "'$COVERAGE_DIR'"$1}'
