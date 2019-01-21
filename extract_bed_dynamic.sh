INPUT=$1
TEMP_DIR=$2
SIZE_FILE=$3

FILENAME=$(basename $INPUT)
FILE_ID=${FILENAME%%.*}

if [ ! -d $TEMP_DIR ]
then
	mkdir -p $TEMP_DIR
fi

echo "INPUT=$INPUT"
echo "TEMP=$TEMP_DIR"
echo "CHROMOSOMES=$SIZE_FILE"
echo "FILE_ID=$FILE_ID"

t1=$(date +%s)
t1_human=$(date)

echo "[STATS] Dividing input table into pieces ["`date`"]"
zcat $INPUT | cut -f 1,2 | awk '{print $0 >> "'$TEMP_DIR'/"$1".table"}'
# read -n 1 -s -r -p "Press any key to continue"

echo "[STATS] Creating single chromosome bed files ["`date`"]"
CHROMOSOMES=()
for chrom in $(cat $SIZE_FILE | cut -f 1)
do
    CHROMOSOMES[${#CHROMOSOMES[@]}]=$chrom
done

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
	if [ -s $TEMP_DIR/$chrom.table ]
	then
	        echo "Calculating bed file for chromosome $chrom = $TEMP_DIR$chrom"
        	python src/cineca/reditools_table_to_bed.py -i $TEMP_DIR/$chrom.table -o $TEMP_DIR/$chrom.bed &
	fi
    done
    wait
    start=$(expr $start + $CHUNK_SIZE)
done
# read -n 1 -s -r -p "Press any key to continue"

t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [BED CHR] [$FILE_ID] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human 1>&2

tmid=$(date +%s)
tmid_human=$(date)

echo "[STATS] Creating complete BED file $TEMP_DIR$FILE_ID.bed ["`date`"]"
rm $TEMP_DIR$FILE_ID".bed"
for chrom in `cat $SIZE_FILE | cut -f 1`
do
	if [ -s $TEMP_DIR$chrom".bed" ]
	then
		cat $TEMP_DIR$chrom".bed" >> $TEMP_DIR$FILE_ID".bed"
	fi
	
	rm $TEMP_DIR$chrom".table"
	rm $TEMP_DIR$chrom".bed"
done

t2=$(date +%s)
t2_human=$(date)
elapsed_time_mid=$(($t2-$tmid))
elapsed_time_mid_human=$(date -d@$elapsed_time_mid -u +%H:%M:%S)
echo "[STATS] [BED GLOBAL] [$FILE_ID] START="$tmid_human" ["$tmid"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time_mid" HUMAN="$elapsed_time_mid_human 1>&2
# read -n 1 -s -r -p "Press any key to continue"

elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)
echo "[STATS] [BED] [$FILE_ID] START="$t1_human" ["$tmid"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human 1>&2

echo -e "$FILE_ID\t$elapsed_time\t$elapsed_time_human" > $TEMP_DIR/$FILE_ID-bed-chronometer.txt

echo "[STATS] Finished creating bed data ["`date`"]"
