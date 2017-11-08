FILENAME=`basename $1`
COVERAGE_DIR=$2
FILE_ID=${FILENAME%.*}

if [ ! -d $COVERAGE_DIR ]
then
	mkdir -p $COVERAGE_DIR
fi

samtools depth $1 | grep -vP "\t0$" | tee $COVERAGE_DIR$FILE_ID".cov" | awk '{print($0) > "'$COVERAGE_DIR'"$1}'
