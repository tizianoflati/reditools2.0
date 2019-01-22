FILENAME=`basename $1`
COVERAGE_DIR=$2
SIZE_FILE=$3

FILE_ID=${FILENAME%.*}

if [ ! -d $COVERAGE_DIR ]
then
	mkdir -p $COVERAGE_DIR
fi

t1=$(date +%s)
t1_human=$(date)

#samtools depth $1 | grep -vP "\t0$" | tee $COVERAGE_DIR$FILE_ID".cov" | awk '{print($0) > "'$COVERAGE_DIR'"$1}'
echo "[STATS] Creating single chromosome coverage files ["`date`"]"
for chrom in `cat $SIZE_FILE | cut -f 1`
do
	echo "Calculating coverage file for chromosome $chrom = $COVERAGE_DIR$chrom"
	
	if [ $(samtools view $1 | cut -f 3 | grep -q $chrom) ]
	then
		samtools depth $1 -r ${chrom#chr} | grep -vP "\t0$" > $COVERAGE_DIR$chrom &
	else
		samtools depth $1 -r $chrom | grep -vP "\t0$" > $COVERAGE_DIR$chrom &
	fi
	
done
wait

t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [COVERAGE CHR] [$FILE_ID] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human 1>&2

tmid=$(date +%s)
tmid_human=$(date)

echo "[STATS] Creating complete file $COVERAGE_DIR$FILE_ID.cov ["`date`"]"
if [ -s $COVERAGE_DIR$FILE_ID".cov" ]
then
	rm $COVERAGE_DIR$FILE_ID".cov"
fi

for chrom in `cat $SIZE_FILE | cut -f 1`
do  
        cat $COVERAGE_DIR$chrom >> $COVERAGE_DIR$FILE_ID".cov"
done

t2=$(date +%s)
t2_human=$(date)
elapsed_time_mid=$(($t2-$tmid))
elapsed_time_mid_human=$(date -d@$elapsed_time_mid -u +%H:%M:%S)
echo "[STATS] [COVERAGE GLOBAL] [$FILE_ID] START="$tmid_human" ["$tmid"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time_mid" HUMAN="$elapsed_time_mid_human 1>&2

elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [COVERAGE] [$FILE_ID] START="$t1_human" ["$tmid"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human 1>&2

echo -e "$FILE_ID\t$elapsed_time\t$elapsed_time_human" > $COVERAGE_DIR/$FILE_ID-coverage-chronometer.txt

echo "[STATS] Finished creating coverage data ["`date`"]"
