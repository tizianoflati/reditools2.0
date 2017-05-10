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
import re

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
        print("[INFO] E[i="+str(key)+"]["+str(len(reads[key]))+"]")
        for read in reads[key]:
#             index = read["index"]
            index = read["alignment_index"]
            
            print("[INFO] \tR:" + str(read["reference"]) + " [pos="+str(read["pos"])+", alignment_index=" + str(index) + ", reference_start="+str(read["object"].reference_start)+" , align_start="+str(read["object"].query_alignment_start) + ", cigar=" + str(read["cigar"])+ ", cigar_list=" + str(read["cigar_list"]) + ", "+ str(len(read["query_qualities"]))+ ", " + str(read["query_qualities"]) + "]")
            print("[INFO] \tQ:" + str(read["sequence"]))
    print("READS[i="+str(i)+"] = " + str(total))

def update_reads(reads):
    if DEBUG:
        print("[INFO] UPDATING READS IN POSITION " + str(i))
    
    for ending_position in reads:
        for read in reads[ending_position]:
#             index = read["index"] - (read["length"] - read["reference_len"])

#             print("ref_len=" +str(read["reference_len"]) + " len="+str(read["length"]) + " diff=" + str(read["reference_len"]-read["length"]));
#             if read["positions"][index] < i:


#             # Original
#             if read["positions"][read["index"]] < i:
#                 read["index"] += 1
#                 read["alignment_index"] += 1
#                 read["reference_index"] += 1

            cigar_list = read["cigar_list"]
            if len(cigar_list) == 0:
                # print("EXCEPTION: CIGAR LIST IS EMPTY")
                continue
            
            if read["pos"] >= i:
                # print("READ POSITION " + str(read["pos"]) + " IS GREATER THAN i=" + str(i))
                continue
            
#             if read["alignment_index"] >= read["length"]:
#                 print("READ FINISHED = " + read["sequence"])
            
            block = cigar_list[0]
            op = block[1]
            
            if op == "S":
                
#                 print(read["object"].reference_start, read)
                
                del cigar_list[0]
#                 if read["sequence"] == "TGGACTTTTCCTGAAATTTATTTTTATGTATGTATATCAAACATTGAATTTCTGTTTTCTTCTTTACTGGAATTGT":
#                     print("[SOFT i="+str(i)+"] S=" + str(block[0])+ " Updating pos from " + str(read["pos"])+ " to " + str(read["pos"] + (block[0])), read["pos"], read)
                
#                 read["pos"] += block[0]
                
#                 if read["reference_index"] >= 0:
#                     read["pos"] += block[0]
#                 else:
#                 read["pos"] += block[0]
                
#                 read["pos"] += block[0]-1
#                 
#                 read["alt"] = read["sequence"][read["alignment_index"]]
#                 read["ref"] = None
#                 #read["qual"] = None
#                 read["qual"] = read["query_qualities"][read["alignment_index"]]
#                 read["alignment_index"] += block[0]
                
                # continue
                
                if not cigar_list:
                    block = None
                else:
                    block = cigar_list[0]
                    op = block[1]
                    
            elif op == "N":
#                 if read["sequence"] == "ATTTTTCTGTTTCTCCCTCAATATCCACCTCATGGAAGTAGATATTCACTAGGTGATATTTTCTAGGCTCTCTTAA":
#                     print("[NNNN i="+str(i)+"] N=" + str(block[0])+ " Updating pos from " + str(read["pos"])+ " to " + str(read["pos"] + (block[0]-1)), read["pos"], read)
                read["pos"] += block[0]
                del cigar_list[0]
                
                read["ref"] = None
                read["alt"] = None
                read["qual"] = None
                
                continue
                
#                 if not cigar_list:
#                     block = None
#                 else:
#                     block = cigar_list[0]
#                     op = block[1]
                    
#             else:
#                 if read["sequence"] == "ATTTTTCTGTTTCTCCCTCAATATCCACCTCATGGAAGTAGATATTCACTAGGTGATATTTTCTAGGCTCTCTTAA":
#                     print("[NORMAL i="+str(i)+"] NORMAL=" + str(block[0])+ " Updating pos from " + str(read["pos"])+ " to " + str(read["pos"] + 1), read["pos"], read)
#                 read["pos"] += 1

            if block is not None and op == "I":
                n = block[0]
                
#                 if read["sequence"] == "GTTAATTTTAGAACATTATCATTCCAAAAAAGCAACTTCATAACATCTAGCAGTCACCTCCTTTCCCATTTCTAGC":
#                     print("[INSERTION i="+str(i)+"] I=" + str(n)+ " Updating alignment_index from " + str(read["alignment_index"]) + " to " + str(read["alignment_index"] + n), read)
                    
#                 read["pos"] += n
                read["alignment_index"] += n
                read["ref"] = None
                read["alt"] = read["sequence"][read["alignment_index"]]
#                 read["alt"] = None
                del cigar_list[0]
                
                if not cigar_list:
                    block = None
                else:
                    block = cigar_list[0]
                    op = block[1]
            
            if block is not None:
                n = block[0]
                
                # D I M N S
                if op == "M":
                    
#                     if read["sequence"] == "GTTAATTTTAGAACATTATCATTCCAAAAAAGCAACTTCATAACATCTAGCAGTCACCTCCTTTCCCATTTCTAGC":
#                         print("[MATCH i="+str(i)+"] M=" + str(n)+ " Updating alignment_index from " + str(read["alignment_index"]) + " to " + str(read["alignment_index"] + 1), read["pos"], read)
                        
                    read["pos"] += 1

                    block[0] -= 1
                    read["reference_index"] += 1
                    read["alignment_index"] += 1
                    
#                     if DEBUG:
#                         print(str(read["reference_index"]), read["reference"][read["reference_index"]], read)
                    
                    read["ref"] = read["reference"][read["reference_index"]]
                    read["alt"] = read["sequence"][read["alignment_index"]]
    
#                     if read["sequence"] == "ATTTTTCTGTTTCTCCCTCAATATCCACCTCATGGAAGTAGATATTCACTAGGTGATATTTTCTAGGCTCTCTTAA":
#                         print("[MATCH i="+str(i)+"]", "pos="+str(read["pos"]), "ref=" + str(read["ref"]), "alt=" + str(read["alt"]), read)
    
                    if block[0] == 0:
                        del cigar_list[0]
                    
                elif op == "D":
#                     if read["sequence"] == "GAAATTTGAAGGTAGAATTGAATACAGATGAACCTCCAATGGTATTCAAGGCTCAGCTGTTTGCGTTGACTGGAGT":
#                         print("[DELETION i="+str(i)+"] D=" + str(n)+ " Updating reference_index from " + str(read["reference_index"])+ " to " + str(read["reference_index"] + n), read["pos"], read)
                    
                    read["reference_index"] += n

                    read["pos"] += n                    
#                     read["alignment_index"] += 1
                    read["ref"] = None
                    # read["ref"] = read["reference"][read["reference_index"]]
                    read["alt"] = None
                    del cigar_list[0]
                
                read["qual"] = read["query_qualities"][read["alignment_index"]]
            
def get_column(reads):
    
    if splice_positions:
        if i in splice_positions[last_chr]:
            if DEBUG: sys.stderr.write("[DEBUG] [SPLICE_SITE] Discarding position ({}, {}) because in splice site\n".format(last_chr, i))
            return None

    if i in omopolymeric_positions:
        if DEBUG:
            sys.stderr.write("[DEBUG] [OMOPOLYMERIC] Discarding position ({}, {}) because omopolymeric\n".format(last_chr, i))
        return None

#     edits = {"T": [], "A": [], "C": [], "G": [], "N": []}    
    edits_no = 0
    edits = []
    ref = None
    
    passed= 0
    qualities = []
    for key in reads:
        for read in reads[key]:
            
#             if DEBUG:
#                 print("GET_COLUMN    i="+str(i) + " READ=" + str(read))
            
            # Filter the reads by positions
            if not filter_base(read):
                
                continue
#             elif read["positions"][read["index"]] != i:
            elif read["pos"] != i:
                if DEBUG:
                    print("[OUT_OF_RANGE] SKIPPING READ i=" + str(i) + " but READ=" + str(read["pos"]))
                continue
            
#             j = read["alignment_index"]
#             if DEBUG:
#                 print("GET_COLUMN_OK i="+str(i) + " ALT="+read["sequence"][j]+" READ=" + str(read))
            
#             ref = read["reference"][read["reference_index"]].upper()
#             if j >= len(read["sequence"]):
#                 print("GET_COLUMN_STRANGE i="+str(i) + " j="+str(j)+" orig="+str(read["alignment_index"])+" READ=" + str(read))
#             alt = read["sequence"][j]

            if read["ref"] == None:
                if DEBUG:
                    print("[INVALID] SKIPPING READ i=" + str(i) + " BECAUSE REF is None", read)
                continue
            if read["alt"] == None:
                if DEBUG:
                    print("[INVALID] SKIPPING READ i=" + str(i) + " BECAUSE ALT is None", read)
                continue
            
            passed += 1

            ref = read["ref"].upper()
            alt = read["alt"].upper()

            edits.append(alt)

            # q = read["query_qualities"][read["alignment_index"]]
            q = read["qual"]
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

    variants = []
#     most_common = None
    ratio = 0.0
    for el in counter.most_common():
        if el[0] == ref: continue
        else:
            variants.append(el[0])
#             most_common = el
            if ratio == 0.0:
                ratio = (float)(el[1]) / (el[1] + ref_count)
        
#     ratio = 0.0
#     if most_common is not None:
#         ratio = (float)(most_common[1]) / (most_common[1] + ref_count)

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
        "variants": variants,
        "frequency": ratio,
        "passed": passed
    }

    # Check that the column passes the filters
    if not filter_column(edits_info): return None

#     if edits_no > 5:
#         print(str(i) + ":" + str(edits_info))
#         raw_input("[ALERT] Press enter to continue...")

    return edits_info;

def filter_read(read):
    
#     if DEBUG:
#         print("[FILTER_READ] F={} QC={} MP={} LEN={} SECOND={} SUPPL={} DUPL={} READ={}".format(read.flag, read.is_qcfail, read.mapping_quality, read.query_length, read.is_secondary, read.is_supplementary, read.is_duplicate, read))
    
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
        if DEBUG: print("[DEBUG] APPLIED FILTER [NOT_PROPER]")
        return False 
 
    return True
    
def filter_base(read):
    
#     pos = read["index"]
    pos = read["alignment_index"]
    
    # Se il carattere e' nelle prime X posizioni della read
    if pos < MIN_BASE_POSITION:
        if DEBUG: print("[DEBUG] APPLIED BASE FILTER [MIN_BASE_POSITION]")
        return False
    
    # Se il carattere e' nelle ultime Y posizioni della read
    if read["length"] - pos < MAX_BASE_POSITION:
        if DEBUG: print("[DEBUG] APPLIED BASE FILTER [MAX_BASE_POSITION]")
        return False
    
    # Se la qualita' e' < Q
    # if read["query_qualities"][read["alignment_index"]] < MIN_BASE_QUALITY:
    if read["qual"] < MIN_BASE_QUALITY:
        if DEBUG: print("[DEBUG] APPLIED BASE FILTER [MIN_BASE_QUALITY]", read["query_qualities"], pos, read["query_qualities"][pos], MIN_BASE_QUALITY, read)
        return False
    
    return True
    
def filter_column(column):
    
    edits = column["edits"]
    
    if column["mean_quality"] < MIN_QUALITY:
        if DEBUG: print("[DEBUG] DISCARDING COLUMN i=" + str(i) + " " + str(column) + " [MIN_MEAN_COLUMN_QUALITY]")
        return False
    
    # Se il numero di caratteri e' < X
    if len(edits) < MIN_COLUMN_LENGTH:
        if DEBUG: print("[DEBUG] DISCARDING COLUMN i=" + str(i) + " " + str(len(edits)) + " [MIN_COLUMN_LENGTH]")
        return False
    
    counter = column["counter"]
    ref = column["ref"]
    
    # (per ogni variazione) se singolarmente il numero delle basi che supportano la variazione e' < X
    for edit in counter:
        if edit != ref and counter[edit] < MIN_EDITS_SINGLE:
            if DEBUG: print("[DEBUG] DISCARDING COLUMN i=" + str(i) + " " + str(counter[edit]) + " [MIN_EDITS_SINGLE]")
            return False
        
    # Se esistono  multipli cambi rispetto al reference
    if len(counter.keys()) > MAX_CHANGES:
        if DEBUG: print("[DEBUG] DISCARDING COLUMN i=" + str(i) + " changes=" + str(len(counter.keys())) + " [MULTIPLE_CHANGES] " + str(column))
        return False
    
    # Se tutte le sostituzioni sono < Y
    if column["edits_no"] < MIN_EDITS_NO:
        if DEBUG: print("[DEBUG] DISCARDING COLUMN i=" + str(i) + " " + str(column["edits_no"]) + " [MIN_EDITS_NO]")
        return False
    
    return True
    
def load_omopolymeric_positions(positions, input_file, region):
    if input_file is None: return
    
    sys.stderr.write("Loading omopolymeric positions from file {}\n".format(input_file))
    
    chromosome = None
    start = None
    end = None
    
    if len(region) >= 1:
        chromosome = region[0]
    if len(region) >= 2:
        start = region[1]
    if len(region) >= 3:
        end = region[2]        
    
    try:
        reader = open(input_file, "r")
        
        for line in reader:
            if line.startswith("#"):
                continue
            
            fields = line.rstrip().split("\t")
            if chromosome is None or fields[0] == chromosome:
                f = int(fields[1])
                t = int(fields[2])
                
                if start is not None: f = max(start, f)
                if end is not None: t = min(t, end)
                
#                 print("POSITION {} {} {} {} {} {}".format(str(fields), chromosome, f, t, start, end))
                
                for i in range(f, t):
                    positions.add(i)
            elif positions:
                break 
            
        reader.close()
    except IOError as e:
        sys.stderr.write("Omopolymeric positions file not found at {}. Error: {}\n".format(input_file, e))
    
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
            sequence = fasta_reader.fetch(chromosome).lower()
            sys.stderr.write("Reference sequence for chromosome {} loaded (len: {})\n".format(chromosome, len(sequence)))
            
            equals = 0
            last = None
            for i, b in enumerate(sequence):
                
#                 if chromosome == "chr18" and i > 190450 and i < 190500:
#                     print(i, b, last, OMOPOLYMERIC_SPAN, sequence[190450:190480])
                
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

def init(samfile, region):
    
    print("Opening bamfile within region=" + str(region))
    
    if len(region) == 0:
        return samfile.fetch()
    
    if len(region) == 1:
        return samfile.fetch(region[0])
    
    else:
        return samfile.fetch(region[0], region[1], region[2])

def within_interval(i, region):
    
    if len(region) <= 1:
        return True
    
    else:
        start = region[1]
        end = region[2]
        return i >= start and i <= end

# -i /marconi_scratch/userexternal/tflati00/test_picardi/reditools_test/SRR1413602.bam
# -o editing18_test -f /marconi_scratch/userinternal/tcastign/test_picardi/hg19.fa -c1,1
# -m20,20 -v1 -q30,30 -e -n0.0 -N0.0 -u -l -p --gzip -H -Y chr18:1-78077248 -F chr18_1_78077248
#
# -f /home/flati/data/reditools/SRR1413602.bam -r /home/flati/data/reditools/hg19.fa -g chr18:14237-14238 -m /home/flati/data/reditools/omopolymeric_positions.txt
#
# -f /home/flati/data/reditools/SRR1413602.bam
# -r /home/flati/data/reditools/hg19.fa
# -g chr18:14237-14238
# -m /home/flati/data/reditools/omopolymeric_positions.txt
if __name__ == '__main__':

    print("START=" + str(datetime.datetime.now()))

    # Options
    MIN_READ_LENGTH = 30 # 100
    SPLICING_SPAN = 5
    OMOPOLYMERIC_SPAN = 5
    MIN_QUALITY = 20 # 20
    MIN_BASE_QUALITY = 30 # 30
    MIN_BASE_POSITION = 6
    MAX_BASE_POSITION = 6 # 76
    MIN_COLUMN_LENGTH = 1 # 10
    MIN_EDITS_SINGLE = 1
    MIN_EDITS_NO = 0
    MAX_CHANGES = 100
    omopolymeric_positions = set()
    
#     MIN_BASE_QUALITY = 20
    MAX_BASE_POSITION = 0
    MIN_BASE_POSITION = 0
    
    # Options parsing
    parser = argparse.ArgumentParser(description='REDItools 2.0')
    parser.add_argument('-f', '--file', help='The bam file to be analyzed')
    parser.add_argument('-o', '--output-file', help='The output statistics file')
    parser.add_argument('-r', '--reference', help='The reference FASTA file')
    parser.add_argument('-g', '--region', help='The region of the bam file to be analyzed')
    parser.add_argument('-m', '--omopolymeric-file', help='The file containing the omopolymeric positions')
    parser.add_argument('-s', '--splicing-file', help='The file containing the splicing sites positions')
    args = parser.parse_args()
    print(args)
    
    bamfile = args.file
    omopolymeric_file = args.omopolymeric_file
    reference_file = args.reference
    output = args.output_file
    splicing_file = args.splicing_file
    
    region = re.split("[:-]", args.region)
    if not region or len(region) == 2 or (len(region) == 3 and region[1] == region[2]):
        sys.stderr.write("[ERROR] Please provide a region of the form chrom:start-end (with end > start). Region provided: {}".format(region))
        exit(1)
    if len(region) >= 2:
        region[1] = int(region[1])
        region[2] = int(region[2])
        
    splice_positions = []
    
    print("Opening BAM file="+bamfile)
    samfile = pysam.AlignmentFile(bamfile, "rb")
    
    load_omopolymeric_positions(omopolymeric_positions, omopolymeric_file, region)
#     if not omopolymeric_positions and omopolymeric_file is not None:
#         omopolymeric_positions = create_omopolymeric_positions(reference_file, omopolymeric_file)
    
    if splicing_file:
        splice_positions = load_splicing_file(splicing_file)
    
    # Constants
    LAST_READ = None
    LOG_INTERVAL = 25000
    
    # Take the time
    tic = datetime.datetime.now()
    
    total = 0
    
    reads = dict()
    reads_list = []
    
    strand = 2
    
    print("Selected region=" + str(region))
    
    outputfile = None
    
    if output is not None:
        outputfile = output
    else:
        prefix = os.path.basename(bamfile)
        if region is not None:
            prefix += "_" + '_'.join([str(x) for x in region])
        outputfile = prefix + "_reditools2_table.gz"
    
    if outputfile.endswith("gz"): writer = gzip.open(outputfile, "w")
    else: writer = open(outputfile, "w")
    
    # Open the iterator
    print("[INFO] Fetching data from bam {}".format(bamfile))
    print("[INFO] Narrowing REDItools to region {}".format(region))
    
    iterator = init(samfile, region)
    
    next_read = next(iterator, LAST_READ)
    if next_read is not None:
        next_pos = next_read.get_reference_positions()
        i = next_pos[0]
        last_chr = next_read.reference_name
        total += 1
        
        read = None
        pos = None
        finished = False
        
        DEBUG_START = -1
        DEBUG_END = -1
        STOP = -1
        
        started = False
        while not finished:
        
            if DEBUG_START > 0 and i >= DEBUG_START: DEBUG = True
            if DEBUG_END > 0 and i >= DEBUG_END: DEBUG = False
            if STOP > 0 and i > STOP: break
        
            if next_read is LAST_READ and len(reads) == 0:
                print("NO MORE READS!")
                finished = True
                break
        
            # Jump if we consumed all the reads
            if len(reads) == 0:
                # print("[INFO] READ SET IS EMPTY. JUMP TO "+str(next_read.get_reference_positions()[0])+"!")
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
                if not filter_read(read): continue
                
                ref_seq = read.get_reference_sequence()
                
    #             if  len(ref_seq) != len(read.query_sequence) or len(pos) != len(read.query_sequence) or len(pos) != len(ref_seq):
    #                 print("=== DETAILS ===")
    #                 print("i="+str(i))
    #                 print("ref_seq="+str(len(ref_seq)))
    #                 print("seq="+str(len(read.query_sequence)))
    #                 print("pos="+str(len(pos)))
    #                 print("qual="+str(len(read.query_qualities)))
    #                 print("index="+str(read.query_alignment_start))
    #                 print(ref_seq)
    #                 print(read.query_sequence)
    #                 print(pos)
    #                 print(read.query_qualities)
                
                item = {
    #                     "index": 0,
                        "pos": read.reference_start - 1,
    #                     "pos": i-1,
                        "alignment_index": read.query_alignment_start - 1,
    #                     "alignment_index": -1,
                        "reference_index": -1,
                        "query_alignment_start": read.query_alignment_start,
                        "object": read,
                        "reference": ref_seq,
                        "reference_len": len(ref_seq),
                        "sequence": read.query_sequence,
                        "positions": pos,
                        "chromosome": read.reference_name,
                        "query_qualities": read.query_qualities,
                        "qualities_len": len(read.query_qualities),
                        "length": read.query_length,
                        "cigar": read.cigarstring
                     }
    
                cigar_list = [[int(c), op] for (c, op) in re.findall('(\d+)(.)', item["cigar"])]
    #             if read.is_reverse:
    #                 cigar_list.reverse()
                item["cigar_list"] = cigar_list
                
    #             if read.query_sequence == "AGGCTCTCTTAATGTAATAAAAGCCATCTATGACAAACCCACAGCCAACATAATACTGAATGGGGAAAAGGTGAAA":
    #                 print(i, read.reference_start, item, read)
                
    #             item["ref"] = item["reference"][item["reference_index"]]
    #             item["alt"] = item["sequence"][item["alignment_index"]]
    #             item["qual"] = item["query_qualities"][item["alignment_index"]]
                
    #             print(item["cigar"])
    #             print(item["cigar_list"])
    #             print(read.get_aligned_pairs())
    #             print("REF START = " + str(read.reference_start))
    #             print("REF POS[0] = " + str(item["positions"][0]))
    #             print("ALIGN START = " + str(item["alignment_index"]))
    #             raw_input("CIGAR STRING PARSED...")
                
    #             if item["cigar"] != "76M":
    #                 item["pairs"] = read.get_aligned_pairs()
                
    #             if read.query_sequence == "CACGGACTTTTCCTGAAATTTATTTTTATGTATGTATATCAAACATTGAATTTCTGTTTTCTTCTTTACTGGAATT" and pos[0] == 14233 and pos[-1] == 14308:
    #                 print("[FILTER_READ] F={} QC={} MP={} LEN={} SECOND={} SUPPL={} DUPL={} PAIRED={} READ={}".format(read.flag, read.is_qcfail, read.mapping_quality, read.query_length, read.is_secondary, read.is_supplementary, read.is_duplicate, read.is_paired, read))
    #                 raw_input("SPECIAL READ...")
                
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
                print("BEFORE UPDATE (i="+str(i)+"):")
                print_reads(reads)
        
            update_reads(reads)
            
            column = get_column(reads)
            
            # Debug purposes
            if DEBUG:
                print("AFTER UPDATE:");
                print_reads(reads)
                raw_input("Press enter to continue...")
            
            # Go the next position
            i += 1
    #         print("Position i"+str(i))
            
            if DEBUG:
                print("[DEBUG] WRITING COLUMN IN POSITION {}: {}".format(i, column is not None))
                print(column)
                print_reads(reads)
            
            if column is not None and within_interval(i, region):
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
                    " ".join([column["ref"] + el for el in column["variants"]]) if column["non_zero"] >= 1 else "-",
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
    