TABLE_DIR=$1
FINAL_TABLE=$2
THREADS=$3

echo "Merging files in $TABLE_DIR using $THREADS threads and writing to output=$FINAL_TABLE"

t1=$(date +%s)
t1_human=$(date)

if [ ! -s $TABLE_DIR/files.txt ]
then
    echo "FILE LIST NOT EXISTING OR EMPTY: "$TABLE_DIR/files.txt
else
    OUTPUT_DIR=`dirname $FINAL_TABLE`
    if [ ! -d $OUTPUT_DIR ]
    then
        mkdir -p $OUTPUT_DIR
    fi
    
    zcat $(cat $TABLE_DIR/files.txt) | bgzip -c -@ $THREADS > $FINAL_TABLE
    echo "Finished creating final table $FINAL_TABLE ["`date`"]"
    
    tabix -s 1 -b 2 -e 2 -c Region $FINAL_TABLE
    echo "Finished creating index file for file $FINAL_TABLE ["`date`"]"    
fi

t2=$(date +%s)
t2_human=$(date)
elapsed_time=$(($t2-$t1))
elapsed_time_human=$(date -d@$elapsed_time -u +%H:%M:%S)

FILE_ID=`basename $TABLE_DIR`

echo "[STATS] [MERGE] [$FILE_ID] START="$t1_human" ["$t1"] END="$t2_human" ["$t2"] ELAPSED="$elapsed_time" HUMAN="$elapsed_time_human 1>&2
echo -e "$FILE_ID\t$elapsed_time\t$elapsed_time_human" > $TABLE_DIR/merge-chronometer.txt

# echo "Starting creating final table "`date`
# 
# i=0; for file in $(ls $TABLE_DIR/*.gz)
# do
#     i=$((i + 1))
#     echo $i". "$file" "`date`; zcat $file >> final_file.txt
# done
# 
# echo "Compressing final_file.txt "`date`
# time /marconi/home/userexternal/tflati00/pigz-2.4/pigz -c final_file.txt > final_file.gz
