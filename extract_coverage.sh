FILENAME=`basename $1`
COVERAGE_DIR=$2
SIZE_FILE=$3

FILE_ID=${FILENAME%.*}

if [ ! -d $COVERAGE_DIR ]
then
	mkdir -p $COVERAGE_DIR
fi

#samtools depth $1 | grep -vP "\t0$" | tee $COVERAGE_DIR$FILE_ID".cov" | awk '{print($0) > "'$COVERAGE_DIR'"$1}'
echo "Creating single chromosome coverage files"
for chrom in `cat $SIZE_FILE | cut -f 1 `
do
	samtools depth $1 -r $chrom | grep -vP "\t0$" > $COVERAGE_DIR$chrom &
done
wait

echo "Creating complete file"
for chrom in `cat $SIZE_FILE | cut -f 1 `
do  
        cat $COVERAGE_DIR$chrom >> $COVERAGE_DIR$FILE_ID".cov"
done
