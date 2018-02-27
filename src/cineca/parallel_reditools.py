#!/usr/bin/env python

import os
import glob
import sys
import re
import time
from mpi4py import MPI
from datetime import datetime
from collections import OrderedDict
import gzip
import reditools
import argparse

ALIGN_CHUNK = 0
STOP_WORKING = 1
IM_FREE = 2
CALCULATE_COVERAGE = 3

STEP = 10000000

TIME_STATS = {}

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

def get_coverage(coverage_file, region = None):

    # Open the file and read i-th section (jump to the next '\n' character)
    n = float(size)
    file_size = os.path.getsize(coverage_file)
    print("[{}] SIZE OF FILE {}: {} bytes".format(rank, coverage_file, file_size))
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
            
            triple = line.rstrip().split("\t")
            
            if region is not None:
                if triple[0] != region[0]: continue
                if len(region) >= 2 and int(triple[1]) < region[1]: continue
                if len(region) >= 2 and int(triple[1]) > region[2]: continue
                
            #if line_no % 10000000 == 0:
            #    print("[{}] [DEBUG] Read {} lines so far".format(rank, line_no))
            cov = int(triple[2])
            coverage_partial += weight_function(cov)

    print("[{}] START={} END={} PARTIAL_COVERAGE={}".format(rank, start, end, coverage_partial))    

    # Reduce
    coverage = None

    coverages = comm.gather(coverage_partial)
    if rank == 0:
        print("COVERAGES:", str(coverages))
        coverage = reduce(lambda x,y: x+y, coverages)
        
    coverage = comm.bcast(coverage, root=0)
        
    # Return the total
    return coverage

def calculate_intervals(total_coverage, coverage_file, region):
    print("[SYSTEM] [{}] Opening coverage file={}".format(rank, coverage_file))
    f = open(coverage_file, "r")

    chr = None
    start = None
    end = None
    C = 0
    max_interval_width = 3000000000 / size

    subintervals = []
    subtotal = total_coverage / size
    print("[SYSTEM] TOTAL={} SUBTOTAL={}".format(total_coverage, subtotal))

    line_no = 0
    with f as lines:
        for line in lines:
            line_no += 1
            if line_no % 1000000 == 0:
                print("[SYSTEM] [{}] Time: {} - {} lines loaded.".format(rank, time.time(), line_no))

            fields = line.rstrip().split("\t")
            
            if region is not None:
                if fields[0] != region[0]: continue
                if len(region) >= 2 and int(fields[1]) < region[1]: continue
                if len(region) >= 3 and int(fields[1]) > region[2]: continue

            # If the interval has become either i) too large or ii) too heavy or iii) spans across two different chromosomes
            if C >= subtotal or (chr is not None and fields[0] != chr) or (end is not None and start is not None and (end-start) > max_interval_width):
                reason = None
                if C >= subtotal: reason = "WEIGHT"
                elif chr is not None and fields[0] != chr: reason = "END_OF_CHROMOSOME"
                elif end is not None and start is not None and (end-start) > max_interval_width: reason = "MAX_WIDTH"

                interval = (chr, start, end, C, end-start, reason)
                print("[SYSTEM] [{}] Time: {} - Discovered new interval={}".format(rank, time.time(), interval))
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
            print("[SYSTEM] [{}] Time: {} - Discovered new interval={}".format(rank, time.time(), interval))
            subintervals.append(interval)

    return subintervals

if __name__ == '__main__':

    # MPI init
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    options = reditools.parse_options()
    options["remove_header"] = True
    
    parser = argparse.ArgumentParser(description='REDItools 2.0')
    parser.add_argument('-G', '--coverage-file', help='The coverage file of the sample to analyze')
    parser.add_argument('-D', '--coverage-dir', help='The coverage directory containing the coverage file of the sample to analyze divided by chromosome')
    parser.add_argument('-t', '--temp-dir', help='The temp directory where to store temporary data for this sample')
    parser.add_argument('-Z', '--chromosome-sizes', help='The file with the chromosome sizes')
    parser.add_argument('-g', '--region', help='The region of the bam file to be analyzed')
    args = parser.parse_known_args()[0]
    
    coverage_file = args.coverage_file
    coverage_dir = args.coverage_dir
    temp_dir = args.temp_dir
    size_file = args.chromosome_sizes
    
    # output = options["output"]
    # format = output.split(".")[-1]
    
    if rank == 0:
        print("[SYSTEM] LAUNCHED PARALLEL REDITOOLS WITH THE FOLLOWING OPTIONS:", options, args)
        
    region = None
    if args.region:
        region = re.split("[:-]", args.region)
        if not region or len(region) == 2 or (len(region) == 3 and region[1] == region[2]):
            sys.stderr.write("[ERROR] Please provide a region of the form chrom:start-end (with end > start). Region provided: {}".format(region))
            exit(1)
        if len(region) >= 2:
            region[1] = int(region[1])
            region[2] = int(region[2])
    
    t1 = time.time()

    print("I am rank #"+str(rank))
    
    if rank == 0:
        print("[0] PRE-COVERAGE TIME " + str(datetime.now().time()))
    
    total_coverage = get_coverage(coverage_file, region)
#     print("TOTAL COVERAGE", str(total_coverage))
    
    if rank == 0:
        now = datetime.now().time()
        elapsed = time.time() - t1
        print("[SYSTEM] [TIME] [MPI] [0] MIDDLE-COVERAGE [now:{}] [elapsed: {}]".format(now, elapsed))
    
    # Collect all the files with the coverage
    files = []
    for file in os.listdir(coverage_dir):
        if region is not None and file != region[0]: continue
        if file.startswith("."): continue
        if file.endswith(".cov"): continue
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
        for i in range(1, min(size, total+1)):
            file = files.pop()
            print("[SYSTEM] [MPI] [0] Sending coverage data "+ str(file) +" to rank " + str(i))
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
            print("[SYSTEM] [TIME] [MPI] [0] COVERAGE RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [#intervals: {}] [{}/{}][{:.2f}%] [Queue:{}]".format(str(who), now, elapsed, len(homeworks), done, total, 100 * float(done)/total, queue))
  
            file = files.pop()
            print("[SYSTEM] [MPI] [0] Sending coverage data "+ str(file) +" to rank " + str(who))
            comm.send(file, dest=who, tag=CALCULATE_COVERAGE)
            queue.add(who)
  
        while len(queue) > 0:
            status = MPI.Status()
            print("[SYSTEM] [MPI] [0] Going to receive data from slaves.")
            subintervals = comm.recv(source=MPI.ANY_SOURCE, tag=IM_FREE, status=status)
            for subinterval in subintervals:
                homeworks.append(subinterval)
  
            done += 1
            who = status.Get_source()
            queue.remove(who)
            now = datetime.now().time()
            elapsed = time.time() - start_intervals
            print("[SYSTEM] [TIME] [MPI] [0] COVERAGE RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [#intervals: {}] [{}/{}][{:.2f}%] [Queue:{}]".format(str(who), now, elapsed, len(homeworks), done, total, 100 * float(done)/total, queue))
  
        now = datetime.now().time()
        elapsed = time.time() - start_intervals
        print("[SYSTEM] [TIME] [MPI] [0] FINISHED CALCULATING INTERVALS [now:{}] [elapsed: {}]".format(now, elapsed))
        
        TIME_STATS["COVERAGE"] = {
                "start": start_intervals,
                "end": now,
                "elapsed": elapsed
            }
        
        done = 0
  
        print("[SYSTEM] [TIME] [MPI] [0] REDItools STARTED. MPI SIZE (PROCS): {} [now: {}]".format(size, datetime.now().time()))
        
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
            
        print("Loading chromosomes' sizes!")
        chromosomes = OrderedDict()
        for line in open(size_file):
            (key, val) = line.split()
            chromosomes[key] = int(val)
        print("Sizes:")
        print(chromosomes)
  
        #homeworks = get_intervals(chromosomes, STEP)
  
        total = len(homeworks)
        print("[SYSTEM] [MPI] [0] HOMEWORKS", total, homeworks)
        #shuffle(homeworks)
  
        start = time.time()
        
        print("[SYSTEM] [TIME] [MPI] [0] REDItools PILEUP START: [now: {}]".format(datetime.now().time()))
        
        queue = set()
        for i in range(1, min(size, total)):
            interval = homeworks.pop()
            print("[SYSTEM] [MPI] [0] Sending data "+ str(interval) +" to rank " + str(i))
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
            print("[SYSTEM] [TIME] [MPI] [0] RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [{}/{}][{:.2f}%] [Queue:{}]".format(str(who), now, elapsed, done, total, 100 * float(done)/total, queue))
  
            interval = homeworks.pop()
            print("[SYSTEM] [MPI] [0] Sending data "+ str(interval) +" to rank " + str(who))
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
            print("[SYSTEM] [TIME] [MPI] [0] RECEIVED IM_FREE SIGNAL FROM RANK {} [now:{}] [elapsed:{}] [{}/{}][{:.2f}%] [Queue:{}]".format(str(who), now, elapsed, done, total, 100 * float(done)/total, queue))
  
        # We have finished processing all the chunks. Let's notify this to slaves
        for i in range(1, size):
            print("[SYSTEM] [MPI] [0] Sending DIE SIGNAL TO RANK " + str(i))
            comm.send(None, dest=i, tag=STOP_WORKING)
  
        #####################################################################
        ######### RECOMBINATION OF SINGLE FILES #############################
        #####################################################################
        t2 = time.time()
        elapsed = t2-t1
        print("[SYSTEM] [TIME] [MPI] [0] WHOLE PARALLEL ANALYSIS FINISHED. CREATING SETUP FOR MERGING PARTIAL FILES - Total elapsed time [{:5.5f}] [{}] [now: {}]".format(elapsed, t2, datetime.now().time()))
        TIME_STATS["COMPUTATION"] = {
                "start": t1,
                "end": t2,
                "elapsed": elapsed
            }
        
        little_files = []
        print("Scanning all files in "+temp_dir+" matching " + ".*")
        for little_file in glob.glob(temp_dir + "/*"):
            print(little_file)
            pieces = re.sub("\..*", "", os.path.basename(little_file)).split("-")
            pieces.insert(0, little_file)
            little_files.append(pieces)

        # Sort the output files            
        print("[SYSTEM] FILES TO MERGE: ", little_files)
        little_files = sorted(little_files, key = lambda x: (x[1], int(x[2])))
        print("[SYSTEM] FILES TO MERGE (SORTED): ", little_files)
        
        smallfiles_list_filename = temp_dir + "files.txt"
        f = open(smallfiles_list_filename, "w")
        for little_file in little_files:
            f.write(f + "\n")
        f.close()
        
        # Open the final output file
#         output_dir = os.path.dirname(output)
#         if not os.path.exists(output_dir):
#             os.makedirs(output_dir)
#         final_file = gzip.open(output, "w")
        
#         final_file.write("\t".join(reditools.get_header()) + "\n")
        
#         total = len(little_files)
#         done = 0
#         for little_file in little_files:
#             print("Writing ", little_file)
#             file = little_file[0]
# 
#             f = gzip.open(file)
#             final_file.write(f.read())
#             f.close()
#             
#             done = done + 1
#             print(file + "\t["+str(done)+"/"+str(total)+" - {:.2%}]".format(done/float(total)))
# 
#         final_file.close()
        
        t2 = time.time()
        print("[SYSTEM] [TIME] [MPI] [0] [END] - WHOLE ANALYSIS FINISHED - Total elapsed time [{:5.5f}] [{}] [now: {}]".format(t2-t1, t2, datetime.now().time()))
        
        print("[STATS] [COVERAGE] START={} END={} ELAPSED={}".format(TIME_STATS["COVERAGE"]["start"], TIME_STATS["COVERAGE"]["end"], TIME_STATS["COVERAGE"]["elapsed"]))
        print("[STATS] [COMPUTATION] START={} END={} ELAPSED={}".format(TIME_STATS["COMPUTATION"]["start"], TIME_STATS["COMPUTATION"]["end"], TIME_STATS["COMPUTATION"]["elapsed"]))
        
    # Slave processes
    if rank > 0:

        while(True):
            # Execute bowtie, view and sort
            status = MPI.Status()
            data = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)

            tag = status.Get_tag()
            if tag == CALCULATE_COVERAGE:
                intervals = calculate_intervals(total_coverage, coverage_dir + data, region)
                comm.send(intervals, dest=0, tag=IM_FREE)
            if tag == ALIGN_CHUNK:

                print("[SYSTEM] [TIME] [MPI] [{}] received data {} from rank 0 [{}]".format(str(rank), str(data), datetime.now().time()))

                # Process it
                time_start = time.time()
                time_s = datetime.now().time()
                print("[SYSTEM] [TIME] [MPI] [{}] REDItools: STARTED [{}]".format(str(rank), time_s))

                # Command: python REDItoolDnaRna_1.04_n.py -i $INPUT -o editing -f hg19.fa -t $THREADS
                # -c 1,1 -m 20,20 -v 0 -q 30,30 -s 2 -g 2 -S -e -n 0.0 -N 0.0 -u -l -H -Y $CHR:$LEFT-$RIGHT -F $CHR_$LEFT_$RIGHT
                # Command REDItools2.0: reditools2.0/src/cineca/reditools.py -f /gss/gss_work/DRES_HAIdA/gtex/SRR1413602/SRR1413602.bam
                #                          -r ../../hg19.fa -g chr18:14237-14238

#                 s = []
# 
#                 if strand == 0:
#                     s = []
#                 elif strand == 1:
#                     s = ['-t', '1', '-i', '2', '-c']
#                 elif strand == 2:
#                     s = ['-t', '2', '-i', '2', '-c']

                id = data[0] + "-" + str(data[1]) + "-" + str(data[2])
                
                options["region"] = [data[0], data[1], data[2]]
                options["output"] = temp_dir + "/" + id + ".gz"
                
#                 command_line = ['time', 'python', 'reditools.py',
#                                 '-f', input,
#                                 '-r', reference,
#                                 '-g', data[0] + ":" + str(data[1]) + "-" + str(data[2]),
#                                 '-m', '/marconi_scratch/userexternal/tflati00/test_picardi/scalability/run/omopolymeric_positions.txt',
#                                 '-o', output + "/" + "-".join([str(rank), data[0], str(data[1]), str(data[2])]) + "." + format]
#                 command_line.extend(s)

#                 print("[MPI] [" + str(rank) + "] COMMAND-LINE:" + ' '.join(command_line))
                print("[MPI] [" + str(rank) + "] COMMAND-LINE:", options)

                reditools.analyze(options)

#                 with open(output+"/logs/error-"+id+".txt","w") as stderr_file:
#                     pid = Popen(command_line, stderr=stderr_file)
#                     pid.wait()

                time_end = time.time()
                print("[SYSTEM] [TIME] [MPI] [{}] REDItools: {} FINISHED [{}][{}] [TOTAL:{:5.2f}]".format(str(rank), str(data), time_s, datetime.now().time(), time_end - time_start))

                print("[SYSTEM] [TIME] [MPI] [{}] sending IM_FREE tag TO RANK 0 [{}]".format(str(rank), datetime.now().time()))
                comm.send(None, dest=0, tag=IM_FREE)
            elif tag == STOP_WORKING:
                print("[SYSTEM] [TIME] [MPI] [{}] received DIE SIGNAL FROM RANK 0 [{}]".format(str(rank), datetime.now().time()))
                break