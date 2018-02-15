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
> pip install -r requirements.txt

Testing
-------------

REDItools2.0 can be tested by issuing the following command

> ./test.sh

Running
-------------

In its most basic form, REDItools 2.0 can be invoked with an input BAM file and an output file:
> python src/cineca/reditools.py -f  $INPUT_BAM_FILE -o table.txt

If you want to restrict the analysis only to a certain region (e.g., only chr1), you can use the **-g** option :
> python src/cineca/reditools.py -f  $INPUT_BAM_FILE -o table.txt -g chr1

Parallel version
-------------

>TODO

Issues
-------------
No known issue are known. For any problem you can write to t.flati@cineca.it.