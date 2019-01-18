# REDItools 2.0

**REDItools 2.0** is the optimized, parallel multi-node version of [<i class="icon-link"></i> REDItools](http://srv00.recas.ba.infn.it/reditools/).


## Installation

### 1. Python setup
---
This guide assumes you have Python installed in your system. If you do not have Python, please read the [official Python webpage](https://www.python.org/).

Make sure you have you preferred Python version loaded. If you have a single Python version already installed in your system you should do nothing. If you have multiple versions, please be sure to point to a given version; in order to do so check your environmental variables (e.g., PATH).

If you are running on a cluster (where usually several versions are available) make sure to load a given Python version. For example (if running on CINECA Marconi super computer) the following command would load Python 2.7.12:
> module load autoload python/2.7.12

Note: REDItools2.0 has been tested with Python 2.7.12. Although there should not be any problem in upgrading, the software comes with no guarantee of being compatible with other versions of Python (e.g., Python >=3).

---
### 2. Cloning / Downloading
---

The first step is to clone this repository (assumes you have *git* installed in your system - see the [Git official page](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) otherwise):
> git clone https://github.com/tflati/reditools2.0.git

(alternatively you can download a ZIP package of REDItools2.0 from [here](https://github.com/tflati/reditools2.0/archive/master.zip) and uncompress the archive).

Move into the project main directory:
> cd reditools2.0

---
### 3. Installing
---

REDItools 2.0 requires a few Python modules to be installed in the environment (e.g., pysam, sortedcontainers, mpi4py, etc.). These can be installed in three ways:

- **System-level**: in this way the dependencies will be installed in your system and all users in your system will see changes. In order to perform this type of installation you need administrator rights.
To install REDItools2.0 in this modality, just run the following command:
> sudo pip install -r requirements.txt

- **User-level**: in this way the dependencies will be installed only for your current user, usually in your home directory. In order to perform this type of installation you need only to be logged as a normal user. Note that this type of installation will install additional software in your local Python directory (usually $HOME/.local/lib/python2.7/site-packages/, but it depends on your operating system and distribution).
This is the recommended modality if you do not care about altering your user environment. Note that altering your user environment might lead to software corruption. For example, assume you have already the *pysam* package installed (version 0.6); since REDItools 2.0 requires a version for *pysam* >= 0.9, the installation would uninstall the existing version of pysam and would install the version 0.9, thus altering the state of your environment. Any existing software which relied on version pysam 0.6 might break and stop working. In conclusion, choose this modality at your own risk.
To install REDItools2.0 in this modality, just run the following command:
> pip install -r requirements.txt --user
 
- **Environment-level**: in this type of installation you create an isolated virtual environment (initially empty) which will contain any new required software, without creating conflicts with any existing environment or requiring any particular right.
This modality will work regardless of the existing packages already installed in your system (both user and system levels) and thus gives the maximum possible freedom to the final-end user.
This is the recommended modality.
The downside of choosing this modality is a potential duplication of code with respect to other existing environments. For example, assume you already have a given version of *sortedcontainers*; by installing REDItools2.0 at environment-level will download and install a *new* copy of *sortedcontainers* into a new isolated environment (ending up with two copies of the same software present in the system, one inside and one outside the virtual environment).
To install REDItools2.0 in this modality, run the following commands:

> virtualenv ENV
source ENV/bin/activate
pip install -r requirements.txt
deactivate

These commands will create a new environment called *ENV* (you can choose any name you like) and will install all dependencies listed in the file *requirements.txt* into it). The commands *activate* and *deactivate* respectively activate (i.e., start/open) and deactivate (i.e., end/close) the virtual environment.

---
### 4. The two versions of REDItools 2.0
---

The software comes with two modalities:
- **Serial version**: in this modality you benefit only from the optimization introduced after the first version. While being significantly faster (with about a 8x factor), you do not exploit the computational power of having multiple cores. On the other hand the setup and launch of REDItools is much easier.
This might be the first modality you might want to give a try when using REDItools2.0 for the first time.

- **Parallel version**: in this modality you benefit both from the serial optimization and from the parallel computation introduced in this brand new version which exploits the existence of multiple cores, also on multiple nodes, making it a perfect tool on High Performance Computing facilities.
Using this modality requires you to perform a little bit more system setup, but it will definitely pay you off.

---
### 5. Testing and running
---
#### 5.1 Serial version

##### Testing
This repo includes test data and a test script for checking that dependencies have been installed properly and the basic REDItools command works.
The serial version of REDItools2.0 can be tested by issuing the following command:

> python src/cineca/reditools.py -f test/SRR2135332.bam -r $REFERENCE -o table.txt -g chr1


##### Running
In its most basic form, REDItools 2.0 can be invoked with an input BAM file, a reference genome and an output file:
> python src/cineca/reditools.py -f \$INPUT_BAM_FILE -r $REFERENCE -o \$OUTPUT_FILE

If you want, you can restrict the analysis only to a certain region (e.g., only chr1), by means of the **-g** option :
> python src/cineca/reditools.py -f  \$INPUT_BAM_FILE -r $REFERENCE -o \$OUTPUT_FILE -g chr1
> 
or a specific interval:
> python src/cineca/reditools.py -f  \$INPUT_BAM_FILE -r $REFERENCE -o \$OUTPUT_FILE -g chr1:1000-2000

For a complete list of options and their usage and meaning, please type:

> python src/cineca/reditools.py -h

---

#### 5.2 Parallel version
The parallel version leverages on the existence of coverage information which reports for each position the number of supporting reads.

We assume you already have 

##### 5.2.1 Producing coverage data
In order to produce such coverage data, execute the script extract_coverage.sh:

> extract_coverage.sh \$FILENAME \$COVERAGE_DIR \$SIZE_FILE

where
$FILENAME is the path of the BAM file to analyze
\$COVERAGE_DIR is the directory that will contain the coverage information
\$SIZE_FILE is the .fai file containing the names of the chromosomes (e.g., hg19.fa.fai).

##### Testing
Testing assumes you have already produced the coverage data (see Section [5.2.1](#521-producing-coverage-data)).
If you can use mpi on your machine (e.g., you are not on a multi-user system and there are no limitations to the jobs you can submit to the system), you can try launching the parallel version of REDItools 2.0 as follows:

> mpirun src/cineca/parallel_reditools.py -f \$SOURCE_BAM_FILE -o \$OUTPUT_DIR/\$SAMPLE_ID/table.gz -r \$REFERENCE -t \$TEMP_DIR -G \$COVERAGE_FILE -D \$COVERAGE_DIR -Z \$SIZE_FILE

This command runs the parallel version of REDItools 2.0 on the input BAM file identified by the variable \$SOURCE_BAM_FILE, saves the output to file \$OUTPUT_DIR/\$SAMPLE_ID/table.gz by using $REFERENCE as the reference genome and using \$TEMP_DIR as temporary directory. The options -G and -D provide the paths of the coverage file and directory, respectively (both created by the *extract_coverage.sh* script).

##### Running

To launch the parallel test on a SLURM-based cluster, just issue the following command:

> sbatch ./parallel_test.sh

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

Issues
-------------
No issues are known so far. For any problem, write to t.flati@cineca.it.
<!--stackedit_data:
eyJoaXN0b3J5IjpbMjA0NjAyNjU3NiwtMjA5NzA0NDIwOCwxMT
U0OTc1MjE0LC05MTM5NDQ4MjNdfQ==
-->