'''
Created on 09 gen 2017

@author: flati
'''

import pysam
import sys
import datetime
from collections import Counter
import gzip
from sortedcontainers import SortedSet
import numpy
import os
import argparse

DEBUG = False

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
            
            print("[INFO] \tR:" + str(read["reference"]) + " [" + str(index) + ", "+ str(len(read["object"].get_reference_positions())) + ", " + str(len(read["query_qualities"]))+ ", " + str(read["object"].get_reference_positions()) + ", " + str(read["query_qualities"]) + "]")
            print("[INFO] \tQ:" + str(read["sequence"]))
    print("READS[i="+str(i)+"] = " + str(total))

def update_reads(pos, reads):
    if DEBUG:
        print("[INFO] UPDATING READS IN POSITION " + str(i))
    
    for ending_position in reads:
        for read in reads[ending_position]:
            index = read["index"] - (read["length"] - read["reference_len"])

#             print("ref_len=" +str(read["reference_len"]) + " len="+str(read["length"]) + " diff=" + str(read["reference_len"]-read["length"]));
            if read["positions"][index] < i:

#             if read["index"] < read["length"]-1:
                read["index"] += 1

def get_column(reads):
    
    if splice_positions:
        if i in splice_positions[last_chr]:
            if DEBUG: sys.stderr.write("[DEBUG] [SPLICE_SITE] Discarding position ({}, {}) because in splice site\n".format(last_chr, i))
            return None

#     edits = {"T": [], "A": [], "C": [], "G": [], "N": []}    
    edits_no = 0
    edits = []
    ref = None
    
    passed= 0
    qualities = []
    for key in reads:
        for read in reads[key]:
            
            # Filter the reads by positions
            if not filter_base(read): continue
            else: passed += 1
            
            j = read["index"]
            
            if DEBUG:
                if j >= len(read["reference"]): 
                    sys.stderr.write("[DEBUG] \t" + str(read["reference"]) + " [" + str(j) + "]")
                    sys.stderr.write("[DEBUG] \t" + str(read["sequence"]))

            length = read["length"]            
            reference_len = read["reference_len"]
            
            ref = read["reference"][j-(length-reference_len)]
            
            ref = ref.upper()
            
            alt = read["sequence"][j]
            edits.append(alt)

            if DEBUG:
                print(read["reference"])
                print(read["sequence"])
                print(j)
                print(ref)
                print(alt)
            
            q = read["query_qualities"][j]
            qualities.append(q)
            
            if alt != ref:
                edits_no += 1
            
    counter = Counter(edits)
    mean_q = 0
    if DEBUG:
        print("Qualities[i="+str(i)+"]="+str(qualities))
        
    if len(qualities) > 0:
        mean_q = numpy.mean(qualities)
    
    # If all the reads are concordant
    #if counter[ref] > 0 and len(counter) == 1:
    #    return None
    
    if len(counter) == 0:
        return None
    
    # [A,C,G,T]
    distribution = [counter['A'], counter['C'], counter['G'], counter['T']]
    ref_count = counter[ref]
    
    non_zero = 0
    for el in counter:
        if el != ref and counter[el] > 0:
            non_zero += 1

    most_common = None
    most_commons = counter.most_common()
    for el in most_commons:
        if el[0] != ref:
            most_common = el
            
    ratio = 0.0
    if most_common is not None:
        ratio = (float)(most_common[1]) / (most_common[1] + ref_count)

#     if passed > 0:
#         print("REF=" + ref)
#         print(passed)
#         print(edits)
#         print(counter)
#         print("MOST FREQUENT EDITS=" + str(counter.most_common()))
#         print("MOST COMMON=" + str(most_common))
#         print(numpy.mean(counter.values()))
#         print(distribution)
#         print(qualities)
#         print(mean_q)
#         print("REF COUNT=" + str(ref_count))
#         print("ALT/REF % = " + str(ratio))
#         raw_input("Press a key:")
    
    edits_info = {
        "edits": edits,
        "distribution": distribution,
        "mean_quality": mean_q,
        "counter": counter,
        "non_zero": non_zero,
        "edits_no": edits_no,
        "ref": ref,
        "most_frequent_edit": counter.most_common(1),
        "frequency": ratio,
        "passed": passed
    }

    # Check that the column passes the filters
    if not filter_column(edits_info): return None

#     if edits_no > 5:
#         print(str(i) + ":" + str(edits_info))
#         raw_input("[ALERT] Press enter to continue...")

    for position in omopolymeric_positions:
        if last_chr == position[0] and i >= position[1] and i <= position[2]:
#             sys.stderr.write("[DEBUG] [OMOPOLYMERIC] Discarding position ({}, {}) because omopolymeric (in region {})\n".format(last_chr, i, position))
            return None

    return edits_info;

def filter_read(read):
    
    # Get the flag of the read
    f = read.flag
    
    # Se la read non e' mappata (FLAG 77 o 141)
    if f == 77 or f == 141:
        if DEBUG: print("[DEBUG] APPLIED FILTER [NOT_MAPPED] f=" + str(f))
        return False 
    
    # Se la read non passa i quality controls (FLAG 512)
    if read.is_qcfail:
        if DEBUG: print("[DEBUG] APPLIED FILTER [QC_FAIL]")
        return False
    
    # Se la read ha un MAPQ < di 30
    if read.mapping_quality < MIN_QUALITY:
        if DEBUG: print("[DEBUG] APPLIED FILTER [MAPQ] " + str(read.mapping_quality) + " MIN="+str(MIN_QUALITY))
        return False
  
    # Se la read ha una lunghezza < XX
    if read.query_length < MIN_READ_LENGTH:
        if DEBUG: print("[DEBUG] APPLIED FILTER [MIN_READ_LENGTH] " + str(read.query_length) + " MIN=" + str(MIN_READ_LENGTH))
        return False

    # Se la read non mappa in modo unico (FLAG 256 o 2048)
    if read.is_secondary or read.is_supplementary:
        if DEBUG: print("[DEBUG] APPLIED FILTER [IS_SECONDARY][IS_SUPPLEMENTARY]")
        return False
    
    # Se la read e' un duplicato di PCR (FLAG 1024)
    if read.is_duplicate:
        if DEBUG: print("[DEBUG] APPLIED FILTER [IS_DUPLICATE]")
        return False

    # Se la read e' paired-end ma non mappa in modo proprio (FLAG diversi da 99/147(+-) o 83/163(-+))
    # 99 = 1+2+32+64 = PAIRED+PROPER_PAIR+MREVERSE+READ1 (+-)
    if read.is_paired and not (f == 99 or f == 147 or f == 83 or f == 163):
        return False 

    return True
    
def filter_base(read):
    
    pos = read["index"]
    
    if DEBUG:
        print(read["query_qualities"][pos], type(read["query_qualities"][pos]), MIN_BASE_QUALITY, read, pos, read["query_qualities"])
    
    # Se il carattere e' nelle prime X posizioni della read
    if pos < MIN_BASE_POSITION:
        if DEBUG: print("[DEBUG] APPLIED BASE FILTER [MIN_BASE_POSITION]")
        return False
    
    # Se il carattere e' nelle ultime Y posizioni della read
    if read["length"] - pos < MAX_BASE_POSITION:
        if DEBUG: print("[DEBUG] APPLIED BASE FILTER [MAX_BASE_POSITION]")
        return False
    
    # Se la qualita' e' < Q
    if read["query_qualities"][pos] < MIN_BASE_QUALITY:
        if DEBUG: print("[DEBUG] APPLIED BASE FILTER [MIN_BASE_QUALITY]", read["query_qualities"], pos, read["query_qualities"][pos], MIN_BASE_QUALITY)
        return False
    
    return True
    
def filter_column(column):
    
    edits = column["edits"]
    
    if column["mean_quality"] < MIN_QUALITY:
        if DEBUG: print("[DEBUG] DISCARDING COLUMN i=" + str(i) + " " + str(column) + " [MIN_MEAN_COLUMN_QUALITY]")
        return False
    
    # TODO: Se il numero di caratteri e' < X
    if len(edits) < MIN_COLUMN_LENGTH:
        if DEBUG: print("[DEBUG] DISCARDING COLUMN i=" + str(i) + " " + str(len(edits)) + " [MIN_COLUMN_LENGTH]")
        return False
    
    counter = column["counter"]
    ref = column["ref"]
    
    # TODO: (per ogni variazione) se singolarmente il numero delle basi che supportano la variazione e' < X
    for edit in counter:
        if edit != ref and counter[edit] < MIN_EDITS_SINGLE:
            if DEBUG: print("[DEBUG] DISCARDING COLUMN i=" + str(i) + " " + str(counter[edit]) + " [MIN_EDITS_SINGLE]")
            return False
        
    # TODO: Se esistono  multipli cambi rispetto al reference
    if len(counter.keys()) > MAX_CHANGES:
        if DEBUG: print("[DEBUG] DISCARDING COLUMN i=" + str(i) + " changes=" + str(len(counter.keys())) + " [MULTIPLE_CHANGES] " + str(column))
        return False
    
    # TODO: Se tutte le sostituzioni sono < Y
    if column["edits_no"] < MIN_EDITS_NO:
        if DEBUG: print("[DEBUG] DISCARDING COLUMN i=" + str(i) + " " + str(column["edits_no"]) + " [MIN_EDITS_NO]")
        return False
    
    return True
    
def load_omopolymeric_positions(positions, input_file):
    if input_file is None: return
    
    sys.stderr.write("Loading omopolymeric positions from file {}\n".format(input_file))
    
    try:
        reader = open(input_file, "r")
        
        for line in reader:
            positions.append(tuple(line.split("\t")))
            
        reader.close()
    except IOError:
        sys.stderr.write("Omopolymeric positions file not found at {}\n".format(input_file))
    
    if not positions:
        sys.stderr.write("Omopolymeric positions file at {} seems to be empty!\n".format(input_file))
    else:
        sys.stderr.write("{} total omopolymeric positions found.\n".format(len(positions)))

def load_chromosome_names(index_file):
    names = []
    
    with open(index_file, "r") as lines:
        for line in lines:
            names.append(line.split("\t")[0])
    
    return names

def load_splicing_file(splicing_file):
    splice_positions = {}
    
    if splicing_file is None: return splice_positions
    
    sys.stderr.write('Loading known splice sites from file {}\n'.format(splicing_file))
    
    if splicing_file.endswith("gz"): f = gzip.open(splicing_file, "r")
    else: f = open(splicing_file, "r")
    
    total = 0
    total_array = {}
    
    for i in f:
            l=(i.strip()).split()
            chr = l[0]
                
            if chr not in splice_positions:
                splice_positions[chr] = SortedSet()
                total_array[chr] = 0           
            
            st,tp,cc = l[4], l[3], int(l[1])
            
            total += SPLICING_SPAN
            total_array[chr] += SPLICING_SPAN
            
            if st=='+' and tp=='D':
                for j in range(SPLICING_SPAN): splice_positions[chr].add(cc+(j+1))
            if st=='+' and tp=='A':
                for j in range(SPLICING_SPAN): splice_positions[chr].add(cc-(j+1))
            if st=='-' and tp=='D':         
                for j in range(SPLICING_SPAN): splice_positions[chr].add(cc-(j+1))
            if st=='-' and tp=='A':
                for j in range(SPLICING_SPAN): splice_positions[chr].add(cc+(j+1))
            
    f.close()
    
    sys.stderr.write('Loaded {} positions from file {}\n'.format(total, splicing_file))
    sys.stderr.write('\tPartial:{}\n'.format(total_array))
    
    return splice_positions

def create_omopolymeric_positions(reference_file, omopolymeric_file):

    tic = datetime.datetime.now()

    sys.stderr.write("Creating omopolymeric positions from reference file {}\n".format(reference_file))
    
    index_file = reference_file + ".fai"
    sys.stderr.write("Loading chromosome names from index file {}\n".format(index_file))
    chromosomes = load_chromosome_names(index_file)
    sys.stderr.write("{} chromosome names found\n".format(str(len(chromosomes))))
    
    positions = []
    
    try:
        # Opening reference fasta file
        sys.stderr.write("Opening reference file {}.\n".format(reference_file))
        fasta_reader = pysam.FastaFile(reference_file)
        sys.stderr.write("Reference file {} opened.\n".format(reference_file))

        for chromosome in chromosomes:
            sys.stderr.write("Loading reference sequence for chromosome {}\n".format(chromosome))
            sequence = fasta_reader.fetch(chromosome)
            sys.stderr.write("Reference sequence for chromosome {} loaded (len: {})\n".format(chromosome, len(sequence)))
                
            equals = 0
            last = None
            for i, b in enumerate(sequence):
                
                if b == last:
                    equals += 1
                else:
                    if equals >= OMOPOLYMERIC_SPAN:
#                         sys.stderr.write("Found a new omopolymeric interval: ({}, {}-{}): {}\n".format(chromosome, i-equals, i, sequence[i-equals:i]))
                        positions.append((chromosome, i-equals, i, equals, last))
                        
                    equals = 1
                     
                last = b
        
        # Closing
        fasta_reader.close()
        sys.stderr.write("Reference file {} closed.\n".format(reference_file))
        
    except ValueError as e:
        sys.stderr.write("Error in reading reference file {}: message={}\n".format(reference_file, e))
    except IOError:
        sys.stderr.write("The reference file {} could not be opened.\n".format(reference_file))

    sys.stderr.write("{} total omopolymeric positions found.\n".format(len(positions)))

    toc = datetime.datetime.now()
    sys.stderr.write("Time to produce all the omopolymeric positions: {}\n".format(toc-tic))

    sys.stderr.write("Writing omopolymeric positions to file: {}.\n".format(omopolymeric_file))
    writer = open(omopolymeric_file, "w")
    writer.write("#" + "\t".join(["Chromomosome", "Start", "End", "Length", "Symbol"]) + "\n")
    for position in positions:
        writer.write("\t".join([str(x) for x in position]) + "\n")
    writer.close()
    sys.stderr.write("Omopolymeric positions written into file: {}.\n".format(omopolymeric_file))

    return positions

# -i /marconi_scratch/userexternal/tflati00/test_picardi/reditools_test/SRR1413602.bam
# -o editing18_test -f /marconi_scratch/userinternal/tcastign/test_picardi/hg19.fa -c1,1 -m20,20 -v1 -q30,30 -e -n0.0 -N0.0 -u -l -p --gzip -H -Y chr18:1-78077248 -F chr18_1_78077248
if __name__ == '__main__':

    print("START=" + str(datetime.datetime.now()))

    # Options
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
    MIN_READ_LENGTH = 30 # 100
    SPLICING_SPAN = 5
    OMOPOLYMERIC_SPAN = 5
    MIN_QUALITY = 30 # 20
    MIN_BASE_QUALITY = 30 # 30
    MIN_BASE_POSITION = 6
    MAX_BASE_POSITION = 6 # 76
    MIN_COLUMN_LENGTH = 1 # 10
    MIN_EDITS_SINGLE = 1
    MIN_EDITS_NO = 0
    MAX_CHANGES = 100
    omopolymeric_positions = []
    
#     MIN_BASE_QUALITY = 20
    MAX_BASE_POSITION = 0
    MIN_BASE_POSITION = 0
    
    # Options parsing
    parser = argparse.ArgumentParser(description='REDItools 2.0')
    parser.add_argument('-f', '--file', help='The bam file to be analyzed')
    parser.add_argument('-r', '--reference', help='The reference FASTA file')
    parser.add_argument('-g', '--region', help='The region of the bam file to be analyzed')
    parser.add_argument('-m', '--omopolymeric-file', help='The file containing the omopolymeric positions')
    parser.add_argument('-s', '--splicing-file', help='The file containing the splicing sites positions')
    args = parser.parse_args()
    print(args)
    
    bamfile = args.file
    omopolymeric_file = args.omopolymeric_file
    reference_file = args.reference
    splicing_file = args.splicing_file
    chromosome_of_interest = args.region
    
    splice_positions = []
    
    print("Opening BAM file="+bamfile)
    samfile = pysam.AlignmentFile(bamfile, "rb")
    
    load_omopolymeric_positions(omopolymeric_positions, omopolymeric_file)
    if not omopolymeric_positions and omopolymeric_file is not None:
        omopolymeric_positions = create_omopolymeric_positions(reference_file, omopolymeric_file)
    
    if splicing_file:
        splice_positions = load_splicing_file(splicing_file)
    
    # Constants
    LAST_READ = None
    LOG_INTERVAL = 100000
    
    # Take the time
    tic = datetime.datetime.now()
    
    total = 0
    
    reads = dict()
    reads_list = []
    
    strand = 2
    
    prefix = os.path.basename(bamfile)
    if chromosome_of_interest is not None:
        prefix += "_" + chromosome_of_interest 
    outputfile = prefix + "_reditools2_table.gz"
    
    if outputfile.endswith("gz"): writer = gzip.open(outputfile, "w")
    else: writer = open(outputfile, "w")
    
    # Open the iterator
    print("[INFO] Fetching data from bam {}".format(bamfile))
    if chromosome_of_interest is not None:
        print("[INFO] Narrowing REDItools to region {}".format(chromosome_of_interest))
        
    iterator = samfile.fetch(chromosome_of_interest)
    
    next_read = next(iterator, LAST_READ)
    next_pos = next_read.get_reference_positions()
    i = next_pos[0]
    last_chr = next_read.reference_name
    total += 1
    
    read = None
    pos = None
    finished = False
    
    started = False
    while not finished:
    
        if i > 14285:
            DEBUG = True
             
        if i > 14287:
            DEBUG = False
            
        if i > 50000:
            break
    
        if next_read is LAST_READ:
                print("NO MORE READS!")
                finished = True
                break
    
        # Jump if we consumed all the reads
        if len(reads) == 0:
            # print("[INFO] READ SET IS EMPTY. JUMP TO "+str(next_read.get_reference_positions()[0])+"!")
            print("setting i=" + str(next_pos[0]))
            i = next_pos[0]
                   
        # Get all the next read(s)
        while next_read is not LAST_READ and (next_pos[0] == i or next_pos[-1] == i):
    
            read = next_read
            pos = next_pos

            next_read = next(iterator, LAST_READ)
            if next_read is not LAST_READ:
                total += 1
                next_pos = next_read.get_reference_positions()
                
                if total % LOG_INTERVAL == 0:
                    print("["+last_chr+"] Total reads loaded: " + str(total) + " ["+str(datetime.datetime.now())+"]")
                
            #print("[INFO] Adding a read to the set=" + str(read.get_reference_positions()))            
            
            # Check that the read passes the filters
#             if not filter_read(read): continue
            
            item = {
                    "index": 0 + (76-len(read.get_reference_sequence())),
                    "object": read,
                    "reference": read.get_reference_sequence(),
                    "reference_len": len(read.get_reference_sequence()),
                    "sequence": read.query_sequence,
                    "positions": pos,
                    "chromosome": read.reference_name,
                    "query_qualities": read.query_qualities,
                    "length": read.query_length
                 }
            
#             print(item)
#             raw_input("Press enter to continue...")
            
#             print(item)
#             if i > 15400000:
#             print(read.reference_name, i)
#             raw_input("Press enter to continue...")
            
            end_position = pos[-1]
            if end_position not in reads:
                reads[end_position] = []
                
            if DEBUG:
                print("Adding item="+str(item))
            reads[end_position].append(item)
            
        # Debug purposes
        if DEBUG:
            print("BEFORE UPDATE:");
            print_reads(reads)
            raw_input("Press enter to continue...")
    
        update_reads(i, reads)
        
        column = get_column(reads)
        
        # Debug purposes
        if DEBUG:
            print("AFTER UPDATE:");
            print_reads(reads)
            raw_input("Press enter to continue...")
        
        # Go the next position
        i += 1
        
        if DEBUG:
            print("[DEBUG] WRITING COLUMN IN POSITION {}: {}".format(i, column is not None))
            print(column)
            print_reads(reads)
        
        if column is not None:
            # head='Region\tPosition\tReference\tStrand\tCoverage-q%i\tMeanQ\tBaseCount[A,C,G,T]\t
            #       AllSubs\tFrequency\t
            #       gCoverage-q%i\tgMeanQ\tgBaseCount[A,C,G,T]\tgAllSubs\tgFrequency\n' %(MQUAL,gMQUAL)
            # cov,bcomp,subs,freq=BaseCount(seq,ref,MINIMUM_EDITS_FREQUENCY,MIN_EDITS_SINGLE)
            # mqua=meanq(qual,len(seq))
            # line='\t'.join([chr,str(pileupcolumn.pos+1),ref,mystrand,str(cov),mqua,str(bcomp),subs,freq]+['-','-','-','-','-'])+'\n'
            writer.write("\t".join([
                last_chr,
                str(i),
                column["ref"],
                str(strand),
                str(column["passed"]),
                "{0:.2f}".format(column["mean_quality"]),
                str(column["distribution"]),
                str(column["ref"] + column["most_frequent_edit"][0][0]) if column["non_zero"] >= 1 else "-",
                "{0:.2f}".format(column["frequency"]),
                "\t".join(['-','-','-','-','-'])
                ]) + "\n")
            writer.flush()
        
        # Remove old reads
        removed = reads.pop(i-1, None)
        
        # When changing chromosome print some statistics
        if read.reference_name != last_chr:
            last_chr = read.reference_name
            
            # if last_chr == "chr2": break
    
            # Take the time
            tac = datetime.datetime.now()
            print("[INFO] REFERENCE NAME=" + last_chr + " (" + str(tac) + ")\t["+delta(tac, tic)+"]")
            tic = tac
    
    samfile.close()
    writer.close()
    
    print("[INFO] TOTAL READS=" + str(total))
    tac = datetime.datetime.now()
    print("[INFO] END=" + str(tac) + "\t["+delta(tac, tic)+"]")
    