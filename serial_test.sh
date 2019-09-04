#!/bin/bash

source ENV/bin/activate

# Serial test #
python src/cineca/reditools.py -f test/SRR2135332.bam -r test/chr21.fa -g chr21 -o serial_table.txt

deactivate
