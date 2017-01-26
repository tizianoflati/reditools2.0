'''
Created on 09 gen 2017

@author: flati
'''

import pysam
import sys
import datetime
import collections

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
    hours, remainder = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
     
    return "%02d:%02d:%02d" % (hours, minutes, seconds)

def print_reads(reads):
    total = 0
    for key in reads:
        total += len(reads[key])
        print("[INFO] E[i="+str(key)+"]")
        for read in reads[key]:
            index = read["index"]
            
            print("[INFO] \t" + str(read["reference"]) + " [" + str(index) + ", "+ str(len(read["object"].get_reference_positions())) + ", " + str(read["object"].get_reference_positions()) + "]")
            print("[INFO] \t" + str(read["sequence"]))
    print("READS[i="+str(i)+"] = " + str(total))

def update_reads(pos, reads):
#     print("[INFO] UPDATING READS IN POSITION " + str(pos))
    
    for ending_position in reads:
        for read in reads[ending_position]:
            if read["positions"][read["index"]] < i:
                read["index"] += 1

def get_edits(reads):
    
    edits_no = 0
#     edits = {"T": [], "A": [], "C": [], "G": [], "N": []}
    edits = []
    edits_info = {"edits": edits}
    
    found = False
    for key in reads:
        for read in reads[key]:
            j = read["index"]
            if j >= len(read["reference"]): 
                print("[DEBUG] \t" + str(read["reference"]) + " [" + str(j) + "]")
                print("[DEBUG] \t" + str(read["sequence"]))
                                
            ref = read["reference"][j]
            if not ref.islower():
                continue
            
            alt = read["sequence"][j]
            edits.append(alt)
            if not found: edits_info["ref"] = ref.upper()
            
            if alt != ref.upper(): edits_no += 1
            
    edits_info["true_edits"] = edits_no
                
#     if edits_no > 5:
#         print(str(i) + ":" + str(edits_info))
#         raw_input("[ALERT] Press enter to continue...")



def filter_read(read):
    # Se la read non e' mappata (FLAG 77)
    if read.is_unmapped: return False
    
    # Se la read ha un MAPQ < di 30
    if read.mapping_quality < 30: return False
    
    # Se la read ha una lunghezza < XX
    if read.query_length < XX: return False
    
    # TODO: Se la read non mappa in modo unico (FLAG 256 o 20148)
#     if read.mapped > 1: return False
    
    # TODO: Se la read e' un duplicato di PCR (FLAG 1024)
    if read.is_duplicate: return False
    
    # TODO: Se la read non passa i quality controls (FLAG 512)
    if read.mapping_quality < MIN_QUALITY: return False
    
    # TODO: Se la read non mappa in modo proprio (FLAG 99/147 o 83/163)
    if read.is_paired and not read.is_proper_pair: return False
    
    return True
    
def filter_base(read, pos):
    
    # TODO: Se il  carattere e' nelle  prime X posizioni della read
    if pos < MIN_BASE_POSITION: return False
    
    # TODO: Se il  carattere e' nelle  ultime Y posizioni della read
    if pos > MAX_BASE_POSITION: return False
    
    # TODO: Se la qualita' e' < Q
    if read.query_qualities[pos] < MIN_BASE_QUALITY: return False
    
    return True
    
def filter_column(column):
    
    edits = column["edits"]
    
    # TODO: Se il numero di caratteri e' < X
    if len(edits) < MIN_COLUMN_LENGTH: return False
    
    counter = collections.Counter(edits)
    ref = column["ref"]
    
    # TODO: (per ogni variazione) se singolarmente il numero delle basi che supportano la variazione e' < X
    for edit in counter:
        if edit != ref and len(counter[edit]) < MIN_EDITS_SINGLE: return False
        
    # TODO: Se esistono  multipli cambi rispetto al reference
    if len(counter.keys()) > 1: return False
    
    # TODO: Se tutte le sostituzioni sono < Y
    if column["true_edits"] < MIN_TRUE_EDITS: return False
    
    return True
    

if __name__ == '__main__':

    print("START=" + str(datetime.datetime.now()))
    
    bamfile = sys.argv[1]
    print("Opening BAM file="+bamfile)
    samfile = pysam.AlignmentFile(bamfile, "rb")
    
    # Options
    XX = 70
    MIN_QUALITY = 100
    MIN_BASE_QUALITY = 200
    MIN_BASE_POSITION = 0
    MAX_BASE_POSITION = 76
    MIN_COLUMN_LENGTH = 10
    MIN_EDITS_SINGLE = 7
    MIN_TRUE_EDITS = 8
    
    # Take the time
    tic = datetime.datetime.now()
    
    total = 0
    last_chr = ""
    
    reads = dict()
    reads_list = []
    
    iterator = samfile.fetch()
    next_read = next(iterator, None)
    next_pos = next_read.get_reference_positions()
    i = next_pos[0]
    total += 1
    
    read = None
    pos = None
    finished = False
    
#     working_chromosome = "4"
            
    started = False
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
                
                if total % 100000 == 0:
                    print("Total reads loaded: " + str(total) + " ["+str(datetime.datetime.now())+"]")
                
            # Check that the read passes the filters
            if not filter_read(read): continue
            
            #print("[INFO] Adding a read to the set=" + str(read.get_reference_positions()))            
            
            item = {
                    "index": 0,
                    "object": read,
                    "reference": read.get_reference_sequence(),
                    "sequence": read.query_sequence,
                    "positions": pos,
                 }
            
            end_position = pos[-1]
            if end_position not in reads:
                reads[end_position] = []
            reads[end_position].append(item)
            
            # Debug purposes
#             print_reads(reads)
#             raw_input("Press enter to continue...")
    
        # Go the next position
        i += 1
        
        # Remove old reads
        removed = reads.pop(i-1, None)
    
        update_reads(i, reads);
        get_edits(reads)

        # When changing chromosome print some statistics
        if read.reference_name != last_chr:
            last_chr = read.reference_name
    
            # Take the time
            tac = datetime.datetime.now()
            print("[INFO] REFERENCE NAME=" + last_chr + " (" + str(tac) + ")\t["+delta(tac, tic)+"]")
            tic = tac
    
    samfile.close()
    
    print("[INFO] TOTAL READS=" + str(total))
    tac = datetime.datetime.now()
    print("[INFO] END=" + str(tac) + "\t["+delta(tac, tic)+"]")
    