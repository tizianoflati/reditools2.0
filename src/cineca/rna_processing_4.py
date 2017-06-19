#!/usr/bin/env python

import os
import glob
import sys
import math
import time
import random
from subprocess import Popen
from mpi4py import MPI
from random import shuffle
from datetime import datetime

ALIGN_CHUNK = 0
STOP_WORKING = 1
IM_FREE = 2
CALCULATE_COVERAGE = 3

STEP = 10000000    

def get_intervals(intervals, num_intervals):
    homeworks = []
    for chromosome in chromosomes.keys():
        print("chromosome:" + chromosome)
        chromosome_length = chromosomes[chromosome]
        print("len:"+ str(chromosome_length))
        step = STEP
        if chromosome == "chrM":
            step = (int)(chromosome_length / 100)
        chromosome_slices = list(range(1, chromosome_length, step)) + [chromosome_length+1]
        print(chromosome_slices)
        print(len(chromosome_slices))
        print("#slices:" + str(len(chromosome_slices)))

        for i in range(0, len(chromosome_slices)-1):
            homeworks.append((chromosome, chromosome_slices[i], chromosome_slices[i+1]-1))
    return homeworks

def weight_function(x):
    # x = math.log(1+x)
    # return 2.748*10**(-3)*x**3 -0.056*x**2 + 0.376*x + 2.093
    return x**3

def get_coverage(coverage_file):

    # Open the file and read i-th section (jump to the next '\n' character)
    n = float(size)
    file_size = os.path.getsize(coverage_file)
    print("[{}] SIZE OF FILE: {}".format(rank, file_size))
    start = int(rank*(file_size/n))
    end = int((rank+1)*(file_size/n))
    print("[{}] [DEBUG] START={} END={}".format(rank, start, end))

    f = open(coverage_file, "r")
    f.seek(start)
    loaded = start
    coverage_partial = 0
    with f as lines:
        line_no = 0
        for line in lines:
            if loaded >= end: continue
            loaded += len(line)

            line_no += 1
            if line_no == 1:
                if not line.startswith("chr"):
                    continue
            #if line_no % 10000000 == 0:
            #    print("[{}] [DEBUG] Read {} lines so far".format(rank, line_no))
            cov = int(line.rstrip().split("\t")[2])
            coverage_partial += weight_function(cov)

    print("[{}] START={} END={} PARTIAL_COVERAGE={}".format(rank, start, end, coverage_partial))    

    # Reduce
    coverage = None

    coverages = comm.gather(coverage_partial)
    if rank == 0:
        print("COVERAGES:", str(coverages))
        coverage = reduce(lambda x,y: x+y, coverages)
        
    # Return the total
    return coverage

def calculate_intervals(total_coverage, coverage_file):
    print("[{}] Opening coverage file={}".format(rank, coverage_file))
    f = open(coverage_file, "r")

    chr = None
    start = None
    end = None
    C = 0
    max_interval_width = 3000000000 / size

    subintervals = []
    subtotal = total_coverage / size
    print("TOTAL={} SUBTOTAL={}".format(total_coverage, subtotal))

    line_no = 0
    with f as lines:
        for line in lines:
            line_no += 1
            if line_no % 1000000 == 0:
                print("[{}] Time: {} - {} lines loaded.".format(rank, time.time(), line_no))

            fields = line.rstrip().split("\t")

            # If the interval has become either i) too large or ii) too heavy or iii) spans across two different chromosomes
            if C >= subtotal or (chr is not None and fields[0] != chr) or (end is not None and start is not None and (end-start) > max_interval_width):
                reason = None
                if C >= subtotal: reason = "WEIGHT"
                elif chr is not None and fields[0] != chr: reason = "END_OF_CHROMOSOME"
                elif end is not None and start is not None and (end-start) > max_interval_width: reason = "MAX_WIDTH"

                interval = (chr, start, end, C, end-start, reason)
                print("[{}] Time: {} - Discovered new interval={}".format(rank, time.time(), interval))
                subintervals.append(interval)
                chr = None
                start = None
                end = None
                C = 0
            if len(fields) < 3: continue

            if chr is None: chr = fields[0]
            if start is None: start = int(fields[1])
            end = int(fields[1])
            # C += math.pow(int(fields[2]), 2)
            C += weight_function(int(fields[2]))

        if C > 0:
            reason = "END_OF_CHROMOSOME"
            interval = (chr, start, end, C, end-start, reason)
            print("[{}] Time: {} - Discovered new interval={}".format(rank, time.time(), interval))
            subintervals.append(interval)

    return subintervals

if __name__ == '__main__':

    # MPI init
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    input=sys.argv[1]
    output=sys.argv[2]
    # strand=sys.argv[3]
    strand = 2
    coverage_file = "coverage_test.txt"

    t1 = time.time()

    print("I am rank #"+str(rank))
    total_coverage = 1324307767144453888
    # total_coverage_cubic_interpolated = 2093347272.03
    # total_coverage = get_coverage(coverage_file)
    # total_coverage_linear = 4909633959
    # total_coverage = 28499712612751
    print("TOTAL COVERAGE", str(total_coverage))
    
    # Collect all the files with the coverage
    files = []
    for file in os.listdir("pieces"):
        if file.startswith("."): continue
        if file == "chrM": continue
        files.append(file)
    files.sort()
    
    if rank == 0:
        print("[0] " + str(len(files)) + " FILES => " + str(files))
    
    '''
    # Assign interval calculation to slaves
    fps = int(len(files) / size)
    if fps == 0: fps = 1
    print("Files per mpi process: " + str(fps))
    subintervals = []
    for i in range(0, size):
        if rank == i:
            from_file = i*fps
            to_file = i*fps+fps if i<size-1 else len(files)
            if from_file > len(files): continue
            if to_file > len(files): continue

            print("[{}] Processing from file {} to file {} = {}".format(rank, from_file, to_file, files[from_file:to_file]))

            for file in files[from_file:to_file]:
                file_intervals = calculate_intervals(total_coverage, "pieces/" + file)
                for interv in file_intervals:
                    subintervals.append(interv)

    # Gather all the intervals calculated from the slaves
    all_subintervals = []
    if rank == 0:
        intervals = None
        all_subintervals = comm.gather(subintervals)
        print("[0] {} total intervals received.".format(len(all_subintervals)))
        homeworks = reduce(lambda x,y: x+y, all_subintervals)
        print("[0] {} total intervals aggregated.".format(len(homeworks)))
        for interval in homeworks:
            print(interval)
    '''

    # Master: dispatches the work to the other slaves
    if rank == 0:
        start_intervals = t1
        print("[0] Start time: {}".format(start_intervals))
        
        done = 0
        total = len(files)
        
	homeworks = []
        queue = set()
        for i in range(1, min(size, total)):
            file = files.pop()
            print("[MPI] [0] Sending coverage data "+ str(file) +" to rank " + str(i))
            comm.send(file, dest=i, tag=CALCULATE_COVERAGE)
            queue.add(i)

        while len(files) > 0:
            status = MPI.Status()
            subintervals = comm.recv(source=MPI.ANY_SOURCE, tag=IM_FREE, status=status)
            for subinterval in subintervals:
                homeworks.append(subinterval)

            done += 1
            who = status.Get_source()
            queue.remove(who)
            now = datetime.now().time()
            elapsed = time.time() - start_intervals
            print("[MPI] [0] COVERAGE RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [#intervals: {}] [{}/{}][{:.2f}%] [Queue:{}]".format(str(who), now, elapsed, len(homeworks), done, total, 100 * float(done)/total, queue))

            file = files.pop()
            print("[MPI] [0] Sending coverage data "+ str(file) +" to rank " + str(who))
            comm.send(file, dest=who, tag=CALCULATE_COVERAGE)
            queue.add(who)

        while len(queue) > 0:
            status = MPI.Status()
            print("[MPI] [0] Going to receive data from slaves.")
            subintervals = comm.recv(source=MPI.ANY_SOURCE, tag=IM_FREE, status=status)
            for subinterval in subintervals:
                homeworks.append(subinterval)

            done += 1
            who = status.Get_source()
            queue.remove(who)
            now = datetime.now().time() 
            elapsed = time.time() - start_intervals
            print("[MPI] [0] COVERAGE RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [#intervals: {}] [{}/{}][{:.2f}%] [Queue:{}]".format(str(who), now, elapsed, len(homeworks), done, total, 100 * float(done)/total, queue))

        print("[MPI] [0] FINISHED CALCULATING INTERVALS")
        done = 0

        print("[MPI] [0] REDItools STARTED")
        print("[MPI] [0] MPI SIZE: " + str(size))
        print(sys.argv)

        if not os.path.exists(output):
            os.makedirs(output)

        print("Loading chromosomes' sizes!")
        chromosomes = {}
        size_file = "/marconi_scratch/userinternal/tcastign/test_picardi/hg19.chrom.sizes";
        with open(size_file) as file:
            for line in file:
                (key, val) = line.split()
                chromosomes[key] = int(val)
        print("Sizes:")
        print(chromosomes)

        #homeworks = get_intervals(chromosomes, STEP)

        total = len(homeworks)
        print("[MPI] [0] HOMEWORKS", total, homeworks)
        #shuffle(homeworks)

        start = time.time()

        queue = set()
        for i in range(1, min(size, total)):
            interval = homeworks.pop()
            print("[MPI] [0] Sending data "+ str(interval) +" to rank " + str(i))
            comm.send(interval, dest=i, tag=ALIGN_CHUNK)
            queue.add(i)

        while len(homeworks) > 0:
            status = MPI.Status()
            comm.recv(source=MPI.ANY_SOURCE, tag=IM_FREE, status=status)
            done += 1
            who = status.Get_source()
            queue.remove(who)
            now = datetime.now().time()
            elapsed = time.time() - start
            print("[MPI] [0] RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [{}/{}][{:.2f}%] [Queue:{}]".format(str(who), now, elapsed, done, total, 100 * float(done)/total, queue))

            interval = homeworks.pop()
            print("[MPI] [0] Sending data "+ str(interval) +" to rank " + str(who))
            comm.send(interval, dest=who, tag=ALIGN_CHUNK)
            queue.add(who)

        while len(queue) > 0:
            status = MPI.Status()
            comm.recv(source=MPI.ANY_SOURCE, tag=IM_FREE, status=status)
            done += 1
            who = status.Get_source()
            queue.remove(who)
            now = datetime.now().time()
            elapsed = time.time() - start
            print("[MPI] [0] RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [{}/{}][{:.2f}%] [Queue:{}]".format(str(who), now, elapsed, done, total, 100 * float(done)/total, queue))

        # We have finished processing all the chunks. Let's notify this to slaves
        for i in range(1, size):
            print("[MPI] [0] Sending DIE SIGNAL TO RANK " + str(i))
            comm.send(None, dest=i, tag=STOP_WORKING)

    # Slave processes
    if rank > 0:

        while(True):
            # Execute bowtie, view and sort
            status = MPI.Status()
            data = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)

            tag = status.Get_tag()
            if tag == CALCULATE_COVERAGE:
                intervals = calculate_intervals(total_coverage, "pieces/" + data)
                comm.send(intervals, dest=0, tag=IM_FREE)
            if tag == ALIGN_CHUNK:

                print("[MPI] [{}] received data {} from rank 0 [{}]".format(str(rank), str(data), datetime.now().time()))

                # Process it
                time_start = time.time()
                time_s = datetime.now().time()
                print("[MPI] [{}] REDItools: STARTED [{}]".format(str(rank), time_s))

                # Command: python REDItoolDnaRna_1.04_n.py -i $INPUT -o editing -f hg19.fa -t $THREADS
                # -c 1,1 -m 20,20 -v 0 -q 30,30 -s 2 -g 2 -S -e -n 0.0 -N 0.0 -u -l -H -Y $CHR:$LEFT-$RIGHT -F $CHR_$LEFT_$RIGHT
                # Command REDItools2.0: reditools2.0/src/cineca/reditools.py -f /gss/gss_work/DRES_HAIdA/gtex/SRR1413602/SRR1413602.bam
                #                          -r ../../hg19.fa -g chr18:14237-14238

                s = []

                if strand == 0:
                    s = []
                elif strand == 1:
                    s = ['-s', '1', '-g', '2', '-S']
                elif strand == 2:
                    s = ['-s', '2', '-g', '2', '-S']

                id = data[0] + "_" + str(data[1]) + "-" + str(data[2])

                command_line = ['time', 'python', '/marconi_scratch/userexternal/tflati00/test_picardi/scalability/run/reditools2.0/src/cineca/reditools.py',
                                '-f', input,
                                '-r', "/marconi_scratch/userexternal/tflati00/test_picardi/hg19.fa",
                                '-g', data[0] + ":" + str(data[1]) + "-" + str(data[2]),
                                '-m', '/marconi_scratch/userexternal/tflati00/test_picardi/scalability/run/omopolymeric_positions.txt',
                                '-o', output + "/" + "-".join([str(rank), data[0], str(data[1]), str(data[2])]) + ".txt"]
                #command_line.extend(s)

                print("[MPI] [" + str(rank) + "] COMMAND-LINE:" + ' '.join(command_line))

                with open(output+"/logs/error-"+id+".txt","w") as stderr_file:
                    pid = Popen(command_line, stderr=stderr_file)
                    pid.wait()

                time_end = time.time()
                print("[MPI] [{}] REDItools: {} FINISHED [{}][{}] [TOTAL:{:5.2f}]".format(str(rank), str(data), time_s, datetime.now().time(), time_end - time_start))

                print("[MPI] [{}] sending IM_FREE tag TO RANK 0 [{}]".format(str(rank), datetime.now().time()))
                comm.send(None, dest=0, tag=IM_FREE)
            elif tag == STOP_WORKING:
                print("[MPI] [{}] received DIE SIGNAL FROM RANK 0 [{}]".format(str(rank), datetime.now().time()))
                break

    comm.Barrier()
    if rank == 0:
        t2 = time.time()
        print("[0] End time: {}".format(t2))
        print("[MPI] [0] WHOLE ANALYSIS FINISHED - Total elapsed time [{:5.5f}]".format(t2-t1))
