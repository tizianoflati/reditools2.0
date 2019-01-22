cd test

if [ ! -s chr21.fa ]
then
	echo "Reference chromosome 21 (Homo Sapiens) not existing. Downloading..."
	wget -O chr21.fa.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz

	echo "Extracting chr21.fa.gz archive"
	gzip -d chr21.fa.gz
fi

if [ ! -s chr21.fa.fai ]
then
	echo "Index .fai not found. Indexing chr21.fa"
	samtools faidx chr21.fa
fi

echo "Test(s) ready!"
