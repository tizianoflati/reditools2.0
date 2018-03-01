TABLE_DIR=$1
FINAL_TABLE=$2
THREADS=$3

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
fi

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
