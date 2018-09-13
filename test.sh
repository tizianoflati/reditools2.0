module load profile/global
module load ig_homo_sapiens/hg19

INPUT_BAM_FILE="test/SRR2135332.bam"
python src/cineca/reditools.py -f $INPUT_BAM_FILE -r $IG_HG19_GENOME/genome.fa -o table.txt -s 1 -g chr1
