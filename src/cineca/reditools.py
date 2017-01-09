'''
Created on 09 gen 2017

@author: flati
'''

import pysam
import sys
import datetime
import time

#'-c', '1,1',
#    INCOV=int(a.split(',')[1])
#    gMINCOV=int(a.split(',')[0])
#'-m', '20,20',
#    MAPQ=int(a.split(',')[1])
#    gMAPQ=int(a.split(',')[0])
#'-v', '1',
#    vnuc=int(a)
#'-q', '30,30',
#    MQUAL=int(a.split(',')[1])
#    gMQUAL=int(a.split(',')[0])
#'-p',
#    conc=1
#'--gzip',
#    gziptab=1
#'-S',
#    corrstr=1
#'-s',
#    getstrand=1
#    if int(a)==1: unchange1,unchange2=1,0
#    elif int(a)==0: unchange1,unchange2=0,0
#    elif int(a)==2: unchange1,unchange2=0,1
#    elif int(a)==12: unchange1,unchange2=1,1
#'-g',
#    if a=='2': useconf=1
#'-e',
#    exh=1
#'-n', '0.0',
#    mmf=float(a)
#'-N', '0.0',
#    gmmf=float(a)
#'-u',
#    mq=1
#'-l',
#    rmsh=1
#'-H',
#    noheader=1

def delta(t2, t1):
    delta = t2 - t1
    hours = 0
    minutes = 0
    seconds = 0
    if hasattr(delta, "hours"):
        hours = delta.hours
    if hasattr(delta, "minutes"):
        minutes = delta.minutes
    if hasattr(delta, "seconds"):
        seconds = delta.seconds
    return "%02d:%02d:%02d" % (hours, minutes, seconds)

def print_reads(reads):
    total = 0
    for key in reads:
        total += len(reads[key])
        #print("[INFO] E[i="+str(key)+"]")
        #for read in reads[key]:
            #pos = read.get_reference_positions()
            #print("[INFO] \t" + str(read.get_reference_sequence()) + " [" + str(pos[0]) + ", " + str(pos[-1]) +"]")
    print("E[i="+str(i)+"] = " + str(total))

if __name__ == '__main__':

    print("START=" + str(datetime.datetime.now()))
    
    bamfile = sys.argv[1]
    print("Opening BAM file="+bamfile)
    samfile = pysam.AlignmentFile(bamfile, "rb")
    
    # Take the time
    tic = datetime.datetime.now()
    
    total = 0
    last_chr = ""
    
    reads = dict()
    
    iterator = samfile.fetch()
    next_read = next(iterator, None)
    next_pos = next_read.get_reference_positions()
    i = next_pos[0]
    total += 1
    
    read = None
    pos = None
    finished = False
    
    #for read in samfile.fetch():
        #total += 1
        #pos = read.get_reference_positions()
        #pos = read.get_reference_positions()
        #pos = read.get_reference_positions()
    
        #if read.reference_name != last_chr:
            #last_chr = read.reference_name
    
            ## Take the time
            #tac = datetime.datetime.now()
            #print("[INFO] REFERENCE NAME=" + last_chr + " (" + str(tac) + ")\t["+delta(tac, tic)+"]")
            #tic = tac
    
    while not finished:
    
        if next_read is None:
                print("NO MORE READS!")
                finished = True
                break
    
        # Jump if we consumed all the reads
        if len(reads) == 0:
            #print("[INFO] READ SET IS EMPTY. JUMP TO "+str(next_read.get_reference_positions()[0])+"!")
            i = next_pos[0]
    
        # Get all the next read(s)
        while next_read is not None and next_pos[0] == i:
    
            read = next_read
            pos = next_pos
    
            next_read = next(iterator, None)
            if next_read is not None:
                total += 1
                next_pos = next_read.get_reference_positions()
    
            #print("[INFO] Adding a read to the set=" + str(read.get_reference_positions()))
            end_position = pos[-1]
            if end_position not in reads:
                reads[end_position] = []
            reads[end_position].append(read)
    
        # Go the next position
        i += 1
    
        # Remove old reads
        removed = reads.pop(i-1, None)
    
        if read.reference_name != last_chr:
            last_chr = read.reference_name
    
            # Take the time
            tac = datetime.datetime.now()
            print("[INFO] REFERENCE NAME=" + last_chr + " (" + str(tac) + ")\t["+delta(tac, tic)+"]")
            tic = tac
    
    samfile.close()
    
    print("[INFO] TOTAL READS=" + str(total))
    print("[INFO] END=" + str(datetime.datetime.now()))
    