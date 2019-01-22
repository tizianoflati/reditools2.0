#!/bin/bash

source ENV/bin/activate

# Serial test #
python src/cineca/reditools.py -f test/SRR2135332.chr21.bam -r test/chr21.fa -o serial_table.txt

deactivate
