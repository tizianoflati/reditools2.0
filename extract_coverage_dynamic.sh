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


echo "[STATS] Creating single chromosome coverage files ["`date`"]"
CHROMOSOMES=()
for chrom in $(cat $SIZE_FILE | cut -f 1)
do
    CHROMOSOMES[${#CHROMOSOMES[@]}]=$chrom
done

###############################
### PER-CHROMOSOME COVERAGE ###
###############################
NUM_CHROMS=$(cat $SIZE_FILE | cut -f 1 | wc -l)
AVAILABLE_CPUS=$(nproc)
CHUNK_SIZE=$(($NUM_CHROMS>$AVAILABLE_CPUS?$AVAILABLE_CPUS:$NUM_CHROMS))
echo "CHROMOSOMES="$NUM_CHROMS
echo "CHUNK SIZE="$CHUNK_SIZE
start=0
while [ $start -lt $NUM_CHROMS ]
do
    echo "NEW BATCH [$(expr $start + 1)-$(expr $start + $CHUNK_SIZE)]"
    for i in $(seq $start $(expr $start + $CHUNK_SIZE - 1))
    do
        if [ $i -ge $NUM_CHROMS ]; then break; fi
        
        chrom=${CHROMOSOMES[$i]}

        echo "Calculating coverage file for chromosome $chrom = $COVERAGE_DIR$chrom"
        
        if [ $(samtools view $1 | cut -f 3 | grep -q $chrom) ]
		then
			samtools depth $1 -r ${chrom#chr} | grep -vP "\t0$" > $COVERAGE_DIR$chrom &
		else
			samtools depth $1 -r $chrom | grep -vP "\t0$" > $COVERAGE_DIR$chrom &
		fi
        
    done
    wait
    start=$(expr $start + $CHUNK_SIZE)
done

t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [COVERAGE CHR] [$FILE_ID] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human 1>&2

tmid=$(date +%s)
tmid_human=$(date)

############################
### SINGLE COVERAGE FILE ###
############################
echo "[STATS] Creating complete file $COVERAGE_DIR$FILE_ID.cov ["`date`"]"
rm $COVERAGE_DIR$FILE_ID".cov"
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
