#!/bin/bash

# Parallel test #
source ENV/bin/activate

SOURCE_BAM_FILE="test/SRR2135332.bam"
REFERENCE="test/chr21.fa"
SIZE_FILE="test/chr21.fa.fai"

NUM_CORES=2
OUTPUT_FILE="test_results/output/parallel_table.txt.gz"
TEMP_DIR="test_results/temp/"
COVERAGE_FILE="test_results/coverage/SRR2135332.chr21.cov"
COVERAGE_DIR="test_results/coverage/"

./extract_coverage.sh $SOURCE_BAM_FILE $COVERAGE_DIR $SIZE_FILE
mpirun -np $NUM_CORES src/cineca/parallel_reditools.py -g "chr21" -f $SOURCE_BAM_FILE -o $OUTPUT_FILE -r $REFERENCE -t $TEMP_DIR -Z $SIZE_FILE -G $COVERAGE_FILE -D $COVERAGE_DIR
./merge.sh $TEMP_DIR $OUTPUT_FILE $NUM_CORES

deactivate
