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

This repo includes test data and a test script for checking that dependencies have been installed properly and the basic REDItools command works.

The software comes with two modalities.

####  4.1 Serial version

In this modality you benefit only from the optimization introduced after the first version. While being significantly faster (with about a 8x factor), you do not exploit the computational power of having multiple cores. On the other hand the setup and launch of REDItools is much easier.
This might be the first modality you might want to give a try when using REDItools2.0 for the first time.

The serial version of REDItools2.0 can be tested by issuing the following command:

> serial_test.sh

or, if you are in a SLURM-based cluster:

> sbatch serial_test_slurm.sh

#### 4.2 Parallel version

In this modality you benefit both from the serial optimization and from the parallel computation introduced in this brand new version which exploits the existence of multiple cores, also on multiple nodes, making it a perfect tool on High Performance Computing facilities.
Using this modality requires you to perform a little bit more system setup, but it will definitely pay you off.

The parallel version leverages on the existence of coverage information which reports for each position the number of supporting reads.

We assume you already have installed and correctly configured the following tools:

- **samtools** (http://www.htslib.org/)
- **htslib** (http://www.htslib.org/)

If you can use *mpi* on your machine (e.g., you are not on a multi-user system and there are no limitations to the jobs you can submit to the system), you can try launching the parallel version of REDItools 2.0 as follows:

> ./parallel_test.sh

If you are running on a SLURM-based cluster, instead, run the following command:

> sbatch ./parallel_test_slurm.sh

This script:
- first defines a bunch of variables which point to input, output and accessory files; then
- launches the production of coverage data; then
- REDItools 2.0 is launched in parallel, by using the specified number of cores; finally
- results are gathered and written into a single table (parameter *-o* provided in the command line)

### 5. Running REDItools 2.0 on your own data

You can now customize the input test scripts to your needs with your input, output and ad-hoc options.


Issues
-------------
No issues are known so far. For any problem, write to t.flati@cineca.it.
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTMxNzI0Mjk0NywyMDQ2MDI2NTc2LC0yMD
k3MDQ0MjA4LDExNTQ5NzUyMTQsLTkxMzk0NDgyM119
-->