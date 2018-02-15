INPUT_BAM_FILE="test/SRR2135332.bam"
python src/cineca/reditools.py -f $INPUT_BAM_FILE -o table.txt -s 1 -g chr1
