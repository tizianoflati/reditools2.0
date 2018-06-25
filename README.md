REDItools 2.0
===================


**REDItools 2.0** is the optimized, parallel multi-node version of [<i class="icon-link"></i> REDItools](http://srv00.recas.ba.infn.it/reditools/).

----------

Installation
-------------
> git clone https://github.com/tflati/reditools2.0
>
> cd reditools2.0
>
> module load autoload python/2.7.12
>
> virtualenv ENV
>
> source ENV/bin/activate
>
> pip install -r requirements.txt
>
>deactivate

or if you plan to use only the serial version:

> git clone https://github.com/tflati/reditools2.0.git
>
> cd reditools2.0
>
> pip install -r requirements.txt

Testing the serial version:
-------------

REDItools2.0 can be tested by issuing the following command

> ./test.sh

Testing the parallel version
-------------

1. Open the parallel test file:
> ./parallel_test.sh 

2. Modify the following variables:

> **BASE_DIR**=\$CINECA_SCRATCH"/reditools/"
> 
> **INPUT_DIR**="/marconi_scratch/userexternal/epicardi/PRJNA231202/SRR1047874/"
> 
> **OUTPUT_DIR**=\$BASE_DIR"/output/"
> 
> **SAMPLE_ID**="SRR1047874"
> 
> **SOURCE_BAM_FILE**=\$INPUT_DIR\$SAMPLE_ID".bam"
> 
> **REFERENCE**=\$BASE_DIR"hg19.fa"
> 
> **OMOPOLYMER_FILE**=\$BASE_DIR"omopolymeric_positions.txt"
> 
> **SIZE_FILE**=\$BASE_DIR"hg19.chrom.sizes"
> 
> **COVERAGE_DIR**=\$BASE_DIR"/cov/"\$SAMPLE_ID"/"
> 
> **COVERAGE_FILE**=\$COVERAGE_DIR\$SAMPLE_ID".cov"
> 
> **TEMP_DIR**=\$BASE_DIR"/temp/"\$SAMPLE_ID"/"
> 
> **strand**=0

3. Launch the parallel test:

> sbatch ./parallel_test.sh

Alternatively, you can launch the raw command:

> mpirun src/cineca/parallel_reditools.py -f \$SOURCE_BAM_FILE -r \$REFERENCE -m \$OMOPOLYMER_FILE -o \$OUTPUT_DIR/\$SAMPLE_ID/table.gz -G \$COVERAGE_FILE -D \$COVERAGE_DIR -t \$TEMP_DIR -Z \$SIZE_FILE \$options | tee $SAMPLE_ID.log

Running the serial version:
-------------

In its most basic form, REDItools 2.0 can be invoked with an input BAM file and an output file:
> python src/cineca/reditools.py -f  \$INPUT_BAM_FILE -o table.txt

If you want to restrict the analysis only to a certain region (e.g., only chr1), you can use the **-g** option :
> python src/cineca/reditools.py -f  \$INPUT_BAM_FILE -o table.txt -g chr1

Issues
-------------
No issues are known so far. For any problem, write to t.flati@cineca.it.
