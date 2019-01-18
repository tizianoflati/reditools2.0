REDItools 2.0
===================


**REDItools 2.0** is the optimized, parallel multi-node version of [<i class="icon-link"></i> REDItools](http://srv00.recas.ba.infn.it/reditools/).

----------

Installation
-------------

Clone this repository:

> git clone https://github.com/tflati/reditools2.0.git
>
> cd reditools2.0

If you are running on a cluster, load a Python module, e.g.:
> module load autoload python/2.7.12

Note: REDItools2.0 has been tested with Python 2.7.12. Although there should not be any problem, the software comes with no guarantee of being compatible with Python 3.

REDItools 2.0 requires a few Python modules to be installed in the environment (e.g., pysam, sortedcontainers, mpi4py, etc.).
It can be installed in three ways:

- **System level**: in this way the dependencies will be installed in your system and all users will see the modifications. In order to perform this type of installation you need administrator rights.

- **User level**: in this way the dependencies will be installed only for your current user, usually in your home directory. In order to perform this type of installation you need only to be logged as a normal user. Note that this type of installation will install additional software in your local Python directory (usually $HOME/.local/lib/python2.7/site-packages/) and is the recommended installation if you do not care about altering your  have any incompatible packages (e.g., different versions of the packages required by REDItools2.0). For example, assume you have already the *pysam* package installed; since REDItools 2.0 requires a version for *pysam* >= 0.9, the installation would fail since the package is already installed in the system but with a different version.
 
- **Environment level**: in this type of installation you create an isolated virtual environment (initially empty) which will contain any new required software, without creating conflicts with your user's existing environment and without any particular right. The downside is a potential duplication of code with regards of other environments. For example, if you already have a given version of *sortedcontainers*, installing REDItools2.0 in this fashion will download and install a new copy of that package inside a new isolated environment (with at least two copies of the same software, one inside and another one outside the virtual environment).

> pip install -r requirements.txt

Otherwise, if you don't have the rights to install the necessary dependencies, you can create an isolated environment:
(Note: if runnning on a cluster

> virtualenv ENV
>
> source ENV/bin/activate
>
> pip install -r requirements.txt
>
> deactivate


Testing the serial version:
-------------

REDItools2.0 can be tested by issuing the following command

> ./test.sh

Testing the parallel version
-------------

By default, the parallel version writes output and temporary directories on the $SCRATCH area, under 'reditools' directory (it will be created if it does not exist). If you wish to modify this settings, open the parallel test file (parallel_test.sh) and modify the following variables as needed:

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

To launch the parallel test on a SLURM-based cluster, just issue the following command:

> sbatch ./parallel_test.sh

In case you are not running on a cluster, you can launch the raw parallel MPI command as follows:

> mpirun src/cineca/parallel_reditools.py -f \$SOURCE_BAM_FILE -r \$REFERENCE -m \$OMOPOLYMER_FILE -o \$OUTPUT_DIR/\$SAMPLE_ID/table.gz -G \$COVERAGE_FILE -D \$COVERAGE_DIR -Z \$SIZE_FILE -t \$TEMP_DIR \$options | tee $SAMPLE_ID.log

Running the serial version:
-------------

In its most basic form, REDItools 2.0 can be invoked with an input BAM file and an output file:
> python src/cineca/reditools.py -f  \$INPUT_BAM_FILE -o table.txt

If you want to restrict the analysis only to a certain region (e.g., only chr1), you can use the **-g** option :
> python src/cineca/reditools.py -f  \$INPUT_BAM_FILE -o table.txt -g chr1

Issues
-------------
No issues are known so far. For any problem, write to t.flati@cineca.it.
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE4MTYyNzQ3NDgsLTkxMzk0NDgyM119
-->