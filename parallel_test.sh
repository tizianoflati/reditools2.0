#!/bin/bash

# Parallel test #
source ENV/bin/activate

NUM_CORES=2
SOURCE_BAM_FILE="test/SRR2135332.chr21.bam"
OUTPUT_FILE="parallel_table.txt"
REFERENCE="test/chr21.fa"
TEMP_DIR="temp/"
SIZE_FILE="test/chr21.fa.fai"
COVERAGE_FILE="coverage/SRR2135332/SRR2135332.cov"
COVERAGE_DIR="coverage/SRR2135332/"

./extract_coverage.sh $SOURCE_BAM_FILE $COVERAGE_DIR $SIZE_FILE
mpirun -np $NUM_CORES src/cineca/parallel_reditools.py -f $SOURCE_BAM_FILE -o $OUTPUT_FILE -r $REFERENCE -t $TEMP_DIR -Z $SIZE_FILE -G $COVERAGE_FILE -D $COVERAGE_DIR
./merge.sh $TEMP_DIR $OUTPUT $NUM_CORES

deactivate
