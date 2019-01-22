#!/usr/bin/env python

'''
Created on 09 gen 2017

@author: flati
'''

import pysam
import sys
import datetime
from collections import defaultdict
import gzip
from sortedcontainers import SortedSet
import os
import argparse
import re
import psutil
import socket
import netifaces

DEBUG = False

def delta(t2, t1):
    delta = t2 - t1
    hours, remainder = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    return "%02d:%02d:%02d" % (hours, minutes, seconds)

def print_reads(reads, i):
    total = 0
    for key in reads:
        total += len(reads[key])
        print("[INFO] E[i="+str(key)+"]["+str(len(reads[key]))+"] strand=" + str(strand))
        for read in reads[key]:
#             index = read["index"]
            index = read["alignment_index"]
            
            print("[INFO] \tR:" + str(read["reference"]) + " [r1="+str(read["object"].is_read1)+", r2="+str(read["object"].is_read2) +", reverse="+str(read["object"].is_reverse) +", pos="+str(read["pos"])+", alignment_index=" + str(index) + ", reference_start="+str(read["object"].reference_start)+" , align_start="+str(read["object"].query_alignment_start) + ", cigar=" + str(read["cigar"])+ ", cigar_list=" + str(read["cigar_list"]) + ", "+ str(len(read["query_qualities"]))+ ", " + str(read["query_qualities"]) + "]")
            print("[INFO] \tQ:" + str(read["sequence"]))
    print("READS[i="+str(i)+"] = " + str(total))

def update_reads(reads, i):
    if DEBUG:
        print("[INFO] UPDATING READS IN POSITION " + str(i))
    
    pos_based_read_dictionary = {}
    
    total = 0
    
    for ending_position in reads:
        for read in reads[ending_position]:

            cigar_list = read["cigar_list"]
            if len(cigar_list) == 0:
                # print("EXCEPTION: CIGAR LIST IS EMPTY")
                continue
            
            if read["pos"] >= i:
                #print("READ POSITION " + str(read["pos"]) + " IS GREATER THAN i=" + str(i))
                continue
            
            total += 1
            
            block = cigar_list[0]
            op = block[1]
            
            if op == "S":
                
                del cigar_list[0]
                
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
                read["qual"] = DEFAULT_BASE_QUALITY
                
                continue
                
            if block is not None and op == "I":
                n = block[0]
                
#                 if read["sequence"] == "GTTAATTTTAGAACATTATCATTCCAAAAAAGCAACTTCATAACATCTAGCAGTCACCTCCTTTCCCATTTCTAGC":
#                     print("[INSERTION i="+str(i)+"] I=" + str(n)+ " Updating alignment_index from " + str(read["alignment_index"]) + " to " + str(read["alignment_index"] + n), read)
                    
                read["alignment_index"] += n
                read["ref"] = None
                read["alt"] = read["sequence"][read["alignment_index"]]
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
                    #if read["reference_index"] >= len(read["reference"]): print("i={} \nSEQ={} \nORG={}".format(read["reference_index"], read["reference"], read["object"].get_reference_sequence()))
                    read["ref"] = read["reference"][read["reference_index"]]
                    read["alt"] = read["sequence"][read["alignment_index"]]
    
#                     if read["sequence"] == "ATTTTTCTGTTTCTCCCTCAATATCCACCTCATGGAAGTAGATATTCACTAGGTGATATTTTCTAGGCTCTCTTAA":
#                         print("[MATCH i="+str(i)+"]", "pos="+str(read["pos"]), "ref=" + str(read["ref"]), "alt=" + str(read["alt"]), read)
    
                    if block[0] == 0:
                        del cigar_list[0]
                    
                elif op == "D":
#                     if read["sequence"] == "GAAATTTGAAGGTAGAATTGAATACAGATGAACCTCCAATGGTATTCAAGGCTCAGCTGTTTGCGTTGACTGGAGT":
#                         print("[DELETION i="+str(i)+"] D=" + str(n)+ " Updating reference_index from " + str(read["reference_index"])+ " to " + str(read["reference_index"] + n), read["pos"], read)
                    
                    #read["reference_index"] += n # MODIFICATO E COMMENTATO IL 26/03/18

                    read["pos"] += n
#                     read["alignment_index"] += 1
                    read["ref"] = None
                    # read["ref"] = read["reference"][read["reference_index"]]
                    read["alt"] = None
                    del cigar_list[0]
                
                if read["query_qualities"] is not None:
                    read["qual"] = read["query_qualities"][read["alignment_index"]]
                
            p = read["pos"]
            if p not in pos_based_read_dictionary: pos_based_read_dictionary[p] = []
            pos_based_read_dictionary[p].append(read)
    
    if DEBUG:
        print("[INFO] READS UPDATED IN POSITION " + str(i) + ":" + str(total))
                
    return pos_based_read_dictionary
            
def get_column(pos_based_read_dictionary, reads, splice_positions, last_chr, omopolymeric_positions, target_positions, i):
    
    if splice_positions:
        if i in splice_positions[last_chr]:
            if VERBOSE:
                sys.stderr.write("[DEBUG] [SPLICE_SITE] Discarding position ({}, {}) because in splice site\n".format(last_chr, i))
            return None

    if omopolymeric_positions:
        if i in omopolymeric_positions[last_chr]:
            if VERBOSE:
                sys.stderr.write("[DEBUG] [OMOPOLYMERIC] Discarding position ({}, {}) because omopolymeric\n".format(last_chr, i))
            return None
        
    if target_positions:
        if (last_chr in target_positions and i not in target_positions[last_chr]) or ("chr"+last_chr in target_positions and i not in target_positions["chr"+last_chr]):
            if VERBOSE:
                sys.stderr.write("[DEBUG] [TARGET POSITIONS] Discarding position ({}, {}) because not in target positions\n".format(last_chr, i))
            return None

#     edits = {"T": [], "A": [], "C": [], "G": [], "N": []}    
    edits_no = 0
    edits = []
    ref = None
    
#     r1r2distribution = Counter()
    r1r2distribution = defaultdict(int)
    
    strand_column = []
    qualities = []
    for key in reads:
        for read in reads[key]:
            
#             if DEBUG:
#                 print("GET_COLUMN  Q_NAME="+ str(read["object"].query_name)+ " READ1=" + str(read["object"].is_read1) + " REVERSE=" + str(read["object"].is_reverse) + " i="+str(i) + " READ=" + str(read))
            
            # Filter the reads by positions
#             if not filter_base(read):
#                 continue

            pos = read["alignment_index"]
    
            # Se il carattere e' nelle prime X posizioni della read
            if pos < MIN_BASE_POSITION:
                if VERBOSE: sys.stderr.write("[DEBUG] APPLIED BASE FILTER [MIN_BASE_POSITION]\n")
                continue
            
            # Se il carattere e' nelle ultime Y posizioni della read
            if read["length"] - pos < MAX_BASE_POSITION:
                if VERBOSE: sys.stderr.write("[DEBUG] APPLIED BASE FILTER [MAX_BASE_POSITION]\n")
                continue
            
            # Se la qualita' e' < Q
            # if read["query_qualities"][read["alignment_index"]] < MIN_BASE_QUALITY:
            if read["qual"] < MIN_BASE_QUALITY:
                if VERBOSE: sys.stderr.write("[DEBUG] APPLIED BASE FILTER [MIN_BASE_QUALITY] {} {} {} {} {}\n".format(str(read["query_qualities"]), pos, str(read["query_qualities"][pos]), MIN_BASE_QUALITY, read))
                continue
            
#             elif read["positions"][read["index"]] != i:
            if read["pos"] != i:
                if DEBUG:
                    print("[OUT_OF_RANGE] SKIPPING READ i=" + str(i) + " but READ=" + str(read["pos"]))
                continue
            
            if DEBUG:
                print("GET_COLUMN  Q_NAME="+ str(read["object"].query_name)+ " READ1=" + str(read["object"].is_read1) + " REVERSE=" + str(read["object"].is_reverse) + " i="+str(i) + " READ=" + str(read))
            
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
            
#             passed += 1
            
#             if passed > 8000:
#                 break

            ref = read["ref"].upper()
            alt = read["alt"].upper()
            
            if DEBUG:
                print("\tBEF={} {}".format(ref, alt))

            if ref == "N" or alt == "N":
                continue
            
            # print(read["pos"], ref, alt, strand, strand == 1, read["object"].is_read1, read["object"].is_read2, read["object"].is_reverse )
            #ref, alt = fix_strand(read, ref, alt)

            if DEBUG:
                print("\tLAT={} {}".format(ref, alt))

            edits.append(alt)

            # q = read["query_qualities"][read["alignment_index"]]
            q = read["qual"]
            qualities.append(q)
            
            strand_column.append(read["strand"])
#             strand_column.append(get_strand(read))
            
            if alt != ref:
                edits_no += 1
                
            r1r2distribution[("R1" if read["object"].is_read1 else "R2") + ("-REV" if read["object"].is_reverse else "")] += 1
    
    if not IS_DNA:
        vstrand = 2
        if strand != 0:
            vstrand = vstand(''.join(strand_column))
            if vstrand == "+": vstrand = 1
            elif vstrand == "-": vstrand = 0
            elif vstrand == "*": vstrand = 2
        
        if vstrand == 0:
            edits = complement_all(edits)
            ref = complement(ref)
        
        if vstrand in [0, 1] and strand_correction:
            edits, strand_column, qualities, qualities_positions = normByStrand(edits, strand_column, qualities, vstrand)

        if DEBUG:
            print(vstrand, ''.join(strand_column))
    else:
        vstrand = "*"
    
    if DEBUG:
        print(r1r2distribution)
#         counter = defaultdict(str)
#         for e in edits: counter[e] += 1
#         print(Counter(edits))
    
#     if i == 62996785:
#         print(edits, strand_column, len(qualities), qualities)
    
    passed = len(edits)
    
#     counter = Counter(edits)
    counter = defaultdict(int)
    for e in edits: counter[e] += 1
    
#     print(Counter(edits), counter)
    
    mean_q = 0
    if DEBUG:
        print("Qualities[i="+str(i)+"]="+str(qualities))
        
    if len(qualities) > 0:
        #mean_q = numpy.mean(qualities)
        mean_q = float(sum(qualities)) / max(len(qualities), 1)
    
    # If all the reads are concordant
    #if counter[ref] > 0 and len(counter) == 1:
    #    return None
    
    if len(counter) == 0:
        if VERBOSE:
            sys.stderr.write("[VERBOSE] [EMPTY] Discarding position ({}, {}) because the associated counter is empty\n".format(last_chr, i))
        return None
    
    # [A,C,G,T]
    distribution = [counter['A'] if 'A' in counter else 0,
                    counter['C'] if 'C' in counter else 0,
                    counter['G'] if 'G' in counter else 0,
                    counter['T'] if 'T' in counter else 0]
    ref_count = counter[ref] if ref in counter else 0
    
    non_zero = 0
    for el in counter:
        if el != ref and counter[el] > 0:
            non_zero += 1

    variants = []
#     most_common = None
    ratio = 0.0
#     most_common = []
#     most_common_value = -1
#     for el in counter:
#         value = counter[el]
#         if value > most_common_value:
#             most_common_value = value
#             most_common = []
#         if value == most_common_value:
#             most_common.append((el, value))

#     for el in Counter(edits).most_common():
    for el in sorted(counter.items(), key=lambda x: x[1], reverse=True):
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
        "passed": passed,
        "strand": vstrand
    }
    
    # Check that the column passes the filters
    if not filter_column(edits_info, i): return None

#     if edits_no > 5:
#         print(str(i) + ":" + str(edits_info))
#         raw_input("[ALERT] Press enter to continue...")

    return edits_info;

def normByStrand(seq_, strand_, squal_, mystrand_):
    
    st='+'
    if mystrand_== 0: st='-'
    seq,strand,qual,squal=[],[],[],''
    for i in range(len(seq_)):
        if strand_[i]==st:
            seq.append(seq_[i])
            strand.append(strand_[i])
            qual.append(squal_[i])
            squal+=chr(squal_[i])
    return seq,strand,qual,squal

# def fix_strand(read, ref, alt):
#     global strand
#     
#     raw_read = read["object"]
#     
#     if (strand == 1 and ((raw_read.is_read1 and raw_read.is_reverse) or (raw_read.is_read2 and not raw_read.is_reverse))) or (strand == 2 and ((raw_read.is_read1 and not raw_read.is_reverse) or (raw_read.is_read2 and raw_read.is_reverse))):
#         return ref, complement(alt)
#     
#     return ref, alt
def get_strand(read):
    global strand
    
    raw_read = read["object"]
    
    if (strand == 1 and ((raw_read.is_read1 and raw_read.is_reverse) or (raw_read.is_read2 and not raw_read.is_reverse))) or (strand == 2 and ((raw_read.is_read1 and not raw_read.is_reverse) or (raw_read.is_read2 and raw_read.is_reverse))):
        return "-"
    
    return "+"            

def filter_read(read):
    
#     if DEBUG:
#         print("[FILTER_READ] F={} QC={} MP={} LEN={} SECOND={} SUPPL={} DUPL={} READ={}".format(read.flag, read.is_qcfail, read.mapping_quality, read.query_length, read.is_secondary, read.is_supplementary, read.is_duplicate, read))
    
    # Get the flag of the read
    f = read.flag

#     if strict_mode:
#         try:
#             NM = read.get_tag("NM")
#             if NM == 0:
# #             print("SKIPPING", MD_value, read.query_sequence, read.reference_start)
#                 return True
#         except KeyError:
#             pass
        
#     if strict_mode:
#         MD = read.get_tag("MD")
# #         print(MD, read.get_reference_sequence(), read.reference_start)
# #         MD = MD.split(":")[1]
#         try:
#             MD_value = int(MD)
# #             print("SKIPPING", MD_value, read.query_sequence, read.reference_start)
#             return True
#         except ValueError:
# #             print("NO MD VALUE")
#             pass
     
    # Se la read non e' mappata (FLAG 77 o 141)
    if f == 77 or f == 141:
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED FILTER [NOT_MAPPED] f={}\n".format(str(f)))
        return False 
     
    # Se la read non passa i quality controls (FLAG 512)
    if read.is_qcfail:
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED FILTER [QC_FAIL]\n")
        return False
     
    # Se la read ha un MAPQ < di 30
    if read.mapping_quality < MIN_QUALITY:
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED FILTER [MAPQ] {} MIN={}\n".format(read.mapping_quality, MIN_QUALITY))
        return False
   
    # Se la read ha una lunghezza < XX
    if read.query_length < MIN_READ_LENGTH:
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED FILTER [MIN_READ_LENGTH] {} MIN={}\n".format(read.query_length, MIN_READ_LENGTH))
        return False
 
    # Se la read non mappa in modo unico (FLAG 256 o 2048)
    if read.is_secondary or read.is_supplementary:
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED FILTER [IS_SECONDARY][IS_SUPPLEMENTARY]\n")
        return False
     
    # Se la read e' un duplicato di PCR (FLAG 1024)
    if read.is_duplicate:
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED FILTER [IS_DUPLICATE]\n")
        return False
 
    # Se la read e' paired-end ma non mappa in modo proprio (FLAG diversi da 99/147(+-) o 83/163(-+))
    # 99 = 1+2+32+64 = PAIRED+PROPER_PAIR+MREVERSE+READ1 (+-)
    if read.is_paired and not (f == 99 or f == 147 or f == 83 or f == 163):
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED FILTER [NOT_PROPER]\n")
        return False
    
    if read.has_tag('SA'):
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED FILTER [CHIMERIC_READ]\n")
        return False
 
    return True
    
def filter_base(read):
    
#     pos = read["index"]
    pos = read["alignment_index"]
    
    # Se il carattere e' nelle prime X posizioni della read
    if pos < MIN_BASE_POSITION:
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED BASE FILTER [MIN_BASE_POSITION]\n")
        return False
    
    # Se il carattere e' nelle ultime Y posizioni della read
    if read["length"] - pos < MAX_BASE_POSITION:
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED BASE FILTER [MAX_BASE_POSITION]\n")
        return False
    
    # Se la qualita' e' < Q
    # if read["query_qualities"][read["alignment_index"]] < MIN_BASE_QUALITY:
    if "qual" not in read:
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED BASE FILTER [QUAL MISSING] {} {}\n".format(pos, read))
        return False
    
    if read["qual"] < MIN_BASE_QUALITY:
        if VERBOSE: sys.stderr.write("[DEBUG] APPLIED BASE FILTER [MIN_BASE_QUALITY] {} {} {} {} {}\n".format(str(read["query_qualities"]), pos, str(read["query_qualities"][pos]), MIN_BASE_QUALITY, read))
        return False
    
    return True
    
def filter_column(column, i):
    
    edits = column["edits"]
    
    if column["mean_quality"] < MIN_QUALITY:
        if VERBOSE: sys.stderr.write("[DEBUG] DISCARDING COLUMN i={} {} [MIN_MEAN_COLUMN_QUALITY]\n".format(i, column))
        return False
    
    # Se il numero di caratteri e' < X
    if len(edits) < MIN_COLUMN_LENGTH:
        if VERBOSE: sys.stderr.write("[DEBUG] DISCARDING COLUMN i={} {} [MIN_COLUMN_LENGTH]\n".format(i, len(edits)))
        return False
    
    counter = column["counter"]
    ref = column["ref"]
    
    # (per ogni variazione) se singolarmente il numero delle basi che supportano la variazione e' < X
    for edit in counter:
        if edit != ref and counter[edit] < MIN_EDITS_SINGLE:
            if VERBOSE: sys.stderr.write("[DEBUG] DISCARDING COLUMN i={} c({})={} [MIN_EDITS_SINGLE] {}\n".format(i, edit, counter[edit], counter))
            return False
        
    # Se esistono  multipli cambi rispetto al reference
    if len(counter.keys()) > MAX_CHANGES:
        if VERBOSE: sys.stderr.write("[DEBUG] DISCARDING COLUMN i={} changes={} [MULTIPLE_CHANGES] {}\n".format(i, len(counter.keys()), column))
        return False
    
    # Se tutte le sostituzioni sono < Y
    if column["edits_no"] < MIN_EDITS_NO:
        if VERBOSE: sys.stderr.write("[DEBUG] DISCARDING COLUMN i={} {} [MIN_EDITS_NO]\n".format(i, column["edits_no"]))
        return False
    
    return True
    
def load_omopolymeric_positions(positions, input_file, region):
    if input_file is None: return
    
    sys.stderr.write("Loading omopolymeric positions from file {}\n".format(input_file))
    
    chromosome = None
    start = None
    end = None
    
    if region is not None:
        if len(region) >= 1:
            chromosome = region[0]
        if len(region) >= 2:
            start = region[1]
        if len(region) >= 3:
            end = region[2]
            
    lines_read = 0
    total = 0
    
    print("Loading omopolymeric positions of {} between {} and {}".format(chromosome, start, end))
    
    try:
        reader = open(input_file, "r")
        
        for line in reader:
            if line.startswith("#"):
                continue
            
            lines_read += 1
            if lines_read % 500000 == 0:
                sys.stderr.write("{} lines read.\n".format(lines_read))
            
            fields = line.rstrip().split("\t")
            if chromosome is None or fields[0] == chromosome:
                chrom = fields[0]
                f = int(fields[1])
                t = int(fields[2])
                
                if start is not None: f = max(start, f)
                if end is not None: t = min(t, end)
                
#                 print("POSITION {} {} {} {} {} {}".format(str(fields), chromosome, f, t, start, end))
                
                if chrom not in positions:
                    positions[chrom] = SortedSet()
                
                for i in range(f, t):
                    positions[chrom].add(i)
                    total += 1
                    
            elif positions:
                break
            
        reader.close()
    except IOError as e:
        sys.stderr.write("[{}] Omopolymeric positions file not found at {}. Error: {}\n".format(region, input_file, e))
    
    sys.stderr.write("[{}] {} total omopolymeric positions found.\n".format(region, total))

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
            chrom = l[0]
                
            if chrom not in splice_positions:
                splice_positions[chrom] = SortedSet()
                total_array[chrom] = 0           
            
            st,tp,cc = l[4], l[3], int(l[1])
            
            total += SPLICING_SPAN
            total_array[chrom] += SPLICING_SPAN
            
            if st=='+' and tp=='D':
                for j in range(SPLICING_SPAN): splice_positions[chrom].add(cc+(j+1))
            if st=='+' and tp=='A':
                for j in range(SPLICING_SPAN): splice_positions[chrom].add(cc-(j+1))
            if st=='-' and tp=='D':         
                for j in range(SPLICING_SPAN): splice_positions[chrom].add(cc-(j+1))
            if st=='-' and tp=='A':
                for j in range(SPLICING_SPAN): splice_positions[chrom].add(cc+(j+1))
            
    f.close()
    
    sys.stderr.write('Loaded {} positions from file {}\n'.format(total, splicing_file))
    sys.stderr.write('\tPartial:{}\n'.format(total_array))
    
    return splice_positions

def create_omopolymeric_positions(reference_file, omopolymeric_file):

    tic = datetime.datetime.now()

    sys.stderr.write("Creating omopolymeric positions (span={}) from reference file {}\n".format(OMOPOLYMERIC_SPAN, reference_file))
    
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

def init(samfile, region):
    
    print("Opening bamfile within region=" + str(region))
    
    if region is None or len(region) == 0:
        return samfile.fetch()
    
    if len(region) == 1:
        try:
            return samfile.fetch(region[0])
        except ValueError:
            return samfile.fetch(region[0].replace("chr", ""))
    
    else:
        try:
            return samfile.fetch(region[0], region[1], region[2])
        except ValueError:
            return samfile.fetch(region[0].replace("chr", ""), region[1], region[2])

def within_interval(i, region):
    
    if region is None or len(region) <= 1:
        return True
    
    else:
        start = region[1]
        end = region[2]
        return i >= start and i <= end

def get_header():
    return ["Region", "Position", "Reference", "Strand", "Coverage-q30", "MeanQ", "BaseCount[A,C,G,T]", "AllSubs", "Frequency", "gCoverage-q30", "gMeanQ", "gBaseCount[A,C,G,T]", "gAllSubs", "gFrequency"]

from collections import Counter
import pickle
def load_target_positions(bed_file, region):
    print("Loading target positions from file {} (region:{})".format(bed_file, region))
    
#     if os.path.exists(bed_file + "save.p"):
#         return pickle.load(open( bed_file + "save.p", "rb" ))
    
    target_positions = {}

    extension = os.path.splitext(bed_file)[1]
    handler = None
    if extension == ".gz":
        handler = gzip.open(bed_file, "r")
    else:
        handler = open(bed_file, "r")
    
    read = 0
    total_positions = 0
    total = Counter()
    with handler as file:
        for line in file:
            read += 1
            fields = line.strip().split("\t")
            chr = fields[0]
            if read % 10000000 == 0: print("[{1}] {0} total lines read. Total positions: {2}".format(read, datetime.datetime.now(), total_positions))
            
            if region != None and chr.replace("chr", "") != region[0].replace("chr", ""): continue
            
            start = int(fields[1])-1
            
            try:
                end = int(fields[2])-1
            except:
                end = start # In case the file has 2 columns only or the third column is not an integer 
            
            intersection_start = max(region[1] if region is not None and len(region)>1 else 0, start)
            intersection_end = min(region[2] if region is not None and len(region)>2 else sys.maxint, end)
            
            
            # If the target region does not intersect the currently analyzed region
            if intersection_end < intersection_start: continue

#             print(line, chr, start, end, intersection_start, intersection_end, total)
            
            # Add target positions
            if chr not in target_positions: target_positions[chr] = SortedSet()
            for i in range(intersection_start, intersection_end+1):
                    
                target_positions[chr].add(i)
                total[chr] += 1
                total_positions += 1
                
    print("### TARGET POSITIONS ###")
    print(total)
    print("TOTAL POSITIONS:", sum(total.values()))
#     pickle.dump(target_positions, open( bed_file + "save.p", "wb" ) )
    
    return target_positions
            
def analyze(options):
    
    global DEBUG
    global activate_debug
    
    print("[SYSTEM]", "PYSAM VERSION", pysam.__version__)
    print("[SYSTEM]", "PYSAM PATH", pysam.__path__)
    
    interface = 'ib0' if 'ib0' in netifaces.interfaces() else netifaces.interfaces()[0]
    hostname = socket.gethostbyaddr(netifaces.ifaddresses(interface)[netifaces.AF_INET][0]['addr'])
    pid = os.getpid()
    hostname_string = hostname[0] + "|" + hostname[2][0] + "|" + str(pid)
    
    bamfile = options["bamfile"]
    region = options["region"]
    reference_file = options["reference"]
    output = options["output"]
    append = options["append"]
    omopolymeric_file = options["omopolymeric_file"]
    splicing_file = options["splicing_file"]
    create_omopolymeric_file = options["create_omopolymeric_file"]
    bed_file = options["bed_file"] if "bed_file" in options else None
    
    LAUNCH_TIME = datetime.datetime.now()
    print("[INFO] ["+str(region)+"] START=" + str(LAUNCH_TIME))

    print("[INFO] Opening BAM file="+bamfile)
    samfile = pysam.AlignmentFile(bamfile, "rb")

    target_positions = {}
    if bed_file is not None:
        target_positions = load_target_positions(bed_file, region)
    
    omopolymeric_positions = {}
    if create_omopolymeric_file is True:
        if omopolymeric_file is not None:
            create_omopolymeric_positions(reference_file, omopolymeric_file)
        else:
            print("[ERROR] You asked to create the omopolymeric file, but you did not specify any output file. Exiting.")
            return
    
    load_omopolymeric_positions(omopolymeric_positions, omopolymeric_file, region)
#     if not omopolymeric_positions and omopolymeric_file is not None:
#         omopolymeric_positions = create_omopolymeric_positions(reference_file, omopolymeric_file)
    
    splice_positions = []
    
    if splicing_file:
        splice_positions = load_splicing_file(splicing_file)
    
    # Constants
    LAST_READ = None
    LOG_INTERVAL = 10000000
    
    # Take the time
    tic = datetime.datetime.now()
    first_tic = tic
    
    total = 0
    
    reads = dict()
    
    outputfile = None
    
    if output is not None:
        outputfile = output
    else:
        prefix = os.path.basename(bamfile)
        if region is not None:
            prefix += "_" + '_'.join([str(x) for x in region])
        outputfile = prefix + "_reditools2_table.gz"
    
    mode = "a" if append else "w"
    
    if outputfile.endswith("gz"): writer = gzip.open(outputfile, mode)
    else: writer = open(outputfile, mode)
    
    if not options["remove_header"]:
        writer.write("\t".join(get_header()) + "\n")
    
    # Open the iterator
    print("[INFO] Fetching data from bam {}".format(bamfile))
    print("[INFO] Narrowing REDItools to region {}".format(region))
    sys.stdout.flush()
    
    reference_reader = None
    if reference_file is not None: reference_reader = pysam.FastaFile(reference_file)
    chr_ref = None
    
    iterator = init(samfile, region)
    
    next_read = next(iterator, LAST_READ)
    if next_read is not None:
#         next_pos = next_read.get_reference_positions()
#         i = next_pos[0]
        i = next_read.reference_start
        
        total += 1
        
        read = None
#         pos = None
        last_chr = None
        finished = False
        
        DEBUG_START = region[1] if region is not None and len(region) > 1 else -1 
        DEBUG_END = region[2] if region is not None and len(region) > 2 else -1
        STOP = -1
        
        while not finished:
        
            if activate_debug and DEBUG_START > 0 and i >= DEBUG_START-1: DEBUG = True
            if activate_debug and DEBUG_END > 0 and i >= DEBUG_END: DEBUG = False
            if STOP > 0 and i > STOP: break
            
#             if i>=46958774:
#                 print(next_read)
#                 print_reads(reads, i)
#                 raw_input()
            
            if (next_read is LAST_READ and len(reads) == 0) or (region is not None and len(region) >= 3 and i > region[2]):
                print("NO MORE READS!")
                finished = True
                break
            
            # Jump if we consumed all the reads
            if len(reads) == 0:
                i = next_read.reference_start
#                 print("[INFO] READ SET IS EMPTY. JUMP TO "+str(next_pos[0])+"!")
#                 if len(next_pos) == 0: i = next_read.reference_start
#                 else: i = next_pos[0]
            
#             print("P1", next_read.query_name, next_pos)
            
            # Get all the next read(s)
            #while next_read is not LAST_READ and (len(next_pos) > 0 and (next_pos[0] == i or next_pos[-1] == i)): # TODO: why or next_pos[-1] == i?
            while next_read is not LAST_READ and next_read.reference_start == i:
            
                read = next_read
#                 pos = next_pos
                
                # When changing chromosome print some statistics
                if read is not LAST_READ and read.reference_name != last_chr:
                    
                    try:
                        chr_ref = reference_reader.fetch(read.reference_name)
                    except KeyError:
                        chr_ref = reference_reader.fetch("chr" + read.reference_name)
                    
                    tac = datetime.datetime.now()
                    print("[INFO] REFERENCE NAME=" + read.reference_name + " (" + str(tac) + ")\t["+delta(tac, tic)+"]")
                    sys.stdout.flush()
                    tic = tac

                last_chr = read.reference_name
                
                next_read = next(iterator, LAST_READ)
                if next_read is not LAST_READ:
                    total += 1
#                     next_pos = next_read.get_reference_positions()
                    
                    if total % LOG_INTERVAL == 0:
                        print("[{}] [{}] [{}] Total reads loaded: {} [{}] [RAM:{}MB]".format(hostname_string, last_chr, region, total, datetime.datetime.now(), psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)))
                        sys.stdout.flush()
                
#                 print("P2", next_read.query_name, next_read.get_reference_positions())
                    
                #print("[INFO] Adding a read to the set=" + str(read.get_reference_positions()))            
                
                # Check that the read passes the filters
                if not filter_read(read): continue
                
#                 ref_seq = read.get_reference_sequence()

                ref_pos = [x[1] for x in read.get_aligned_pairs() if x[0] is not None and x[1] is not None]
                ref_seq = ''.join([chr_ref[x] for x in ref_pos]).upper()
                
#                 if ref_seq != read.get_reference_sequence().upper():
#                     print("MY_REF={} \nPY_REF={} \nREAD_NAME={} \nPOSITIONS={} \nREAD={}\n--------------------------".format(ref_seq, read.get_reference_sequence().upper(), read.query_name, read.get_reference_positions(), read.query_sequence))
                    
#                 raw_input()
                
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
                
                t = "*"
                
                if not IS_DNA:
                    if read.is_read1:
                        if strand == 1:
                            if read.is_reverse: t='-'
                            else: t='+'
                        else:
                            if read.is_reverse: t='+'
                            else: t='-'
                    elif read.is_read2:
                        if strand == 2:
                            if read.is_reverse: t='-'
                            else: t='+'
                        else:
                            if read.is_reverse: t='+'
                            else: t='-'
                    else: # for single ends
                        if strand == 1:
                            if read.is_reverse: t='-'
                            else: t='+'
                        else:
                            if read.is_reverse: t='+'
                            else: t='-'                        
                
                qualities = read.query_qualities
                if qualities is None: qualities = [DEFAULT_BASE_QUALITY for x in range(0, len(ref_seq))]
                
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
#                         "positions": pos,
                        "chromosome": read.reference_name,
                        "query_qualities": qualities,
                        "qualities_len": len(qualities),
                        "length": read.query_length,
                        "cigar": read.cigarstring,
                        "strand": t
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
                
                end_position = read.reference_end #pos[-1]
                if end_position not in reads:
                    reads[end_position] = []
                    
                if DEBUG:
                    print("Adding item="+str(item))
                reads[end_position].append(item)
                
            # Debug purposes
#             if DEBUG:
#                 print("BEFORE UPDATE (i="+str(i)+"):")
#                 print_reads(reads, i)
        
            pos_based_read_dictionary = update_reads(reads, i)
            
            column = get_column(pos_based_read_dictionary, reads, splice_positions, last_chr, omopolymeric_positions, target_positions, i)
            
            # Debug purposes
            if DEBUG:
#                 print("AFTER UPDATE:");
#                 print_reads(reads, i)
                raw_input("Press enter to continue...")
            
            # Go the next position
            i += 1
    #         print("Position i"+str(i))
            
#             if DEBUG:
#                 print("[DEBUG] WRITING COLUMN IN POSITION {}: {}".format(i, column is not None))
#                 print(column)
#                 print_reads(reads, i)
            
            if column is not None and within_interval(i, region) and not (strict_mode and column["non_zero"] == 0):
                # head='Region\tPosition\tReference\tStrand\tCoverage-q%i\tMeanQ\tBaseCount[A,C,G,T]\t
                #       AllSubs\tFrequency\t
                #       gCoverage-q%i\tgMeanQ\tgBaseCount[A,C,G,T]\tgAllSubs\tgFrequency\n' %(MQUAL,gMQUAL)
                # cov,bcomp,subs,freq=BaseCount(seq,ref,MINIMUM_EDITS_FREQUENCY,MIN_EDITS_SINGLE)
                # mqua=meanq(qual,len(seq))
                # line='\t'.join([chr,str(pileupcolumn.pos+1),ref,mystrand,str(cov),mqua,str(bcomp),subs,freq]+['-','-','-','-','-'])+'\n'
                # [A,C,G,T]
                
                writer.write("\t".join([
                    last_chr,
                    str(i),
                    column["ref"],
                    str(column["strand"]),
                    str(column["passed"]),
                    "{0:.2f}".format(column["mean_quality"]),
                    str(column["distribution"]),
                    " ".join([column["ref"] + el for el in column["variants"]]) if column["non_zero"] >= 1 else "-",
                    "{0:.2f}".format(column["frequency"]),
                    "\t".join(['-','-','-','-','-'])
                    ]) + "\n")
#                 if column["passed"] >= 1000: print("WRITTEN LINE {} {} {} {} {}".format(last_chr, str(i), column["ref"], column["strand"], column["passed"]))
                # writer.flush()
            elif VERBOSE:
                sys.stderr.write("[VERBOSE] [NOPRINT] Not printing position ({}, {}) WITHIN_INTERVAL={} STRICT_MODE={} COLUMN={}\n".format(last_chr, i, within_interval(i, region), strict_mode, column))
            
            # Remove old reads
            reads.pop(i-1, None)
    
    if reference_reader is not None: reference_reader.close()
    samfile.close()
    writer.close()
    
    tac = datetime.datetime.now()
    print("[INFO] ["+hostname_string+"] ["+str(region)+"] " + str(total) + " total reads read")
    print("[INFO] ["+hostname_string+"] ["+str(region)+"] END=" + str(tac) + "\t["+delta(tac, tic)+"]")
    print("[INFO] ["+hostname_string+"] ["+str(region).ljust(50)+"] FINAL END=" + str(tac) + " START="+ str(first_tic) + "\t"+ str(region) +"\t[TOTAL COMPUTATION="+delta(tac, first_tic)+"] [LAUNCH TIME:"+str(LAUNCH_TIME)+"] [TOTAL RUN="+delta(tac, LAUNCH_TIME)+"] [READS="+str(total)+"]")

complement_map = {"A":"T", "T":"A", "C":"G", "G":"C"}
def complement(b):
    return complement_map[b]

def complement_all(sequence):
    return ''.join([complement_map[l] for l in sequence])

def prop(tot,va):
    try: av=float(va)/tot
    except: av=0.0
    return av

def vstand(strand): # strand='+-+-+-++++++-+++'
    
    vv=[(strand.count('+'),'+'),(strand.count('-'),'-'),(strand.count('*'),'*')]
    if vv[0][0]==0 and vv[1][0]==0: return '*'
    if use_strand_confidence: #flag che indica se usare il criterio 2, altrimenti usa il criterio 1
        totvv=sum([x[0] for x in vv[:2]])
        if prop(totvv,vv[0][0])>=strand_confidence_value: return '+' # strand_confidence_value e' il valore soglia, compreso tra 0 e 1, default 0.7
        if prop(totvv,vv[1][0])>=strand_confidence_value: return '-'
        return '*'
    else:
        if vv[0][0]==vv[1][0] and vv[2][0]==0: return '+'
        return max(vv)[1]

def parse_options():
    
    # Options parsing
    parser = argparse.ArgumentParser(description='REDItools 2.0')
    parser.add_argument('-f', '--file', help='The bam file to be analyzed')
    parser.add_argument('-o', '--output-file', help='The output statistics file')
    parser.add_argument('-S', '--strict', default=False, action='store_true', help='Activate strict mode: only sites with edits will be included in the output')
    parser.add_argument('-s', '--strand', type=int, default=0, help='Strand: this can be 0 (unstranded), 1 (secondstrand oriented) or 2 (firststrand oriented)')
    parser.add_argument('-a', '--append-file', action='store_true', help='Appends results to file (and creates if not existing)')
    parser.add_argument('-r', '--reference', help='The reference FASTA file')
    parser.add_argument('-g', '--region', help='The region of the bam file to be analyzed')
    parser.add_argument('-m', '--omopolymeric-file', help='The file containing the omopolymeric positions')
    parser.add_argument('-c', '--create-omopolymeric-file', default=False, help='Whether to create the omopolymeric span', action='store_true')
    parser.add_argument('-os', '--omopolymeric-span', type=int, default=5, help='The omopolymeric span')
    parser.add_argument('-sf', '--splicing-file', help='The file containing the splicing sites positions')
    parser.add_argument('-ss', '--splicing-span', type=int, default=4, help='The splicing span')
    parser.add_argument('-mrl', '--min-read-length', type=int, default=30, help='The minimum read length. Reads whose length is below this value will be discarded.')
    parser.add_argument('-q', '--min-read-quality', type=int, default=20, help='The minimum read quality. Reads whose mapping quality is below this value will be discarded.')
    parser.add_argument('-bq', '--min-base-quality', type=int, default=30, help='The minimum base quality. Bases whose quality is below this value will not be included in the analysis.')
    parser.add_argument('-mbp', '--min-base-position', type=int, default=0, help='The minimum base position. Bases which reside in a previous position (in the read) will not be included in the analysis.')
    parser.add_argument('-Mbp', '--max-base-position', type=int, default=0, help='The maximum base position. Bases which reside in a further position (in the read) will not be included in the analysis.')
    parser.add_argument('-l', '--min-column-length', type=int, default=1, help='The minimum length of editing column (per position). Positions whose columns have length below this value will not be included in the analysis.')
    parser.add_argument('-men', '--min-edits-per-nucleotide', type=int, default=1, help='The minimum number of editing for events each nucleotide (per position). Positions whose columns have bases with less than min-edits-per-base edits will not be included in the analysis.')
    parser.add_argument('-me', '--min-edits', type=int, default=0, help='The minimum number of editing events (per position). Positions whose columns have bases with less than \'min-edits-per-base edits\' will not be included in the analysis.')    
    parser.add_argument('-Men', '--max-editing-nucleotides', type=int, default=100, help='The maximum number of editing nucleotides, from 0 to 4 (per position). Positions whose columns have more than \'max-editing-nucleotides\' will not be included in the analysis.')
    parser.add_argument('-d', '--debug', default=False, help='REDItools is run in DEBUG mode.', action='store_true')
    parser.add_argument('-T', '--strand-confidence', default=1, help='Strand inference type 1:maxValue 2:useConfidence [1]; maxValue: the most prominent strand count will be used; useConfidence: strand is assigned if over a prefixed frequency confidence (-TV option)')
    parser.add_argument('-C', '--strand-correction', default=False, help='Strand correction. Once the strand has been inferred, only bases according to this strand will be selected.', action='store_true')
    parser.add_argument('-Tv', '--strand-confidence-value', type=float, default=0.7, help='Strand confidence [0.70]')    
    parser.add_argument('-V', '--verbose', default=False, help='Verbose information in stderr', action='store_true')
    parser.add_argument('-H', '--remove-header', default=False, help='Do not include header in output file', action='store_true')
    parser.add_argument('-N', '--dna', default=False, help='Run REDItools 2.0 on DNA-Seq data', action='store_true')
    parser.add_argument('-B', '--bed_file', help='Path of BED file containing target regions')
    
    args = parser.parse_known_args()[0]
#     print(args)
    
    global activate_debug
    activate_debug = args.debug
    
    global VERBOSE
    VERBOSE = args.verbose
    
    bamfile = args.file
    if bamfile is None:
        print("[ERROR] An input bam file is mandatory. Please, provide one (-f|--file)")
        exit(1)
        
    omopolymeric_file = args.omopolymeric_file
    global OMOPOLYMERIC_SPAN
    OMOPOLYMERIC_SPAN = args.omopolymeric_span
    create_omopolymeric_file = args.create_omopolymeric_file
    
    reference_file = args.reference
    if reference_file is None:
        print("[ERROR] An input reference file is mandatory. Please, provide one (-r|--reference)")
        exit(1)
    
    output = args.output_file
    append = args.append_file
    
    global strict_mode
    strict_mode = args.strict
    
    global strand
    strand = args.strand
    
    global strand_correction
    strand_correction = args.strand_correction
    
    global use_strand_confidence
    use_strand_confidence = bool(args.strand_confidence)
    
    global strand_confidence_value
    strand_confidence_value = float(args.strand_confidence_value)
    
    splicing_file = args.splicing_file
    global SPLICING_SPAN
    SPLICING_SPAN = args.splicing_span
    
    global MIN_READ_LENGTH 
    MIN_READ_LENGTH = args.min_read_length
    
    global MIN_QUALITY
    MIN_QUALITY = args.min_read_quality
    
    global MIN_BASE_QUALITY
    MIN_BASE_QUALITY = args.min_base_quality
    
    global DEFAULT_BASE_QUALITY
    DEFAULT_BASE_QUALITY = 30 
    
    global MIN_BASE_POSITION
    MIN_BASE_POSITION = args.min_base_position
    
    global MAX_BASE_POSITION
    MAX_BASE_POSITION = args.max_base_position
    
    global MIN_COLUMN_LENGTH
    MIN_COLUMN_LENGTH = args.min_column_length
    
    global MIN_EDITS_SINGLE
    MIN_EDITS_SINGLE = args.min_edits_per_nucleotide
    
    global MIN_EDITS_NO
    MIN_EDITS_NO = args.min_edits
    
    global MAX_CHANGES
    MAX_CHANGES = args.max_editing_nucleotides
    
    global IS_DNA
    IS_DNA = args.dna
    
    bed_file = args.bed_file
    
    if IS_DNA and bed_file is None:
        print("[ERROR] When analyzing DNA-Seq files it is mandatory to provide a BED file containing the positions of target regions (-B|--bed_file)")
        exit(1)
    
    region = None
    
    if args.region:
        region = re.split("[:-]", args.region)
        if not region or len(region) == 2 or (len(region) == 3 and region[1] == region[2]):
            sys.stderr.write("[ERROR] Please provide a region of the form chrom:start-end (with end > start). Region provided: {}".format(region))
            exit(1)
        if len(region) >= 2:
            region[1] = int(region[1])
            region[2] = int(region[2])
    
    options = {
        "bamfile": bamfile,
        "region": region,
        "reference": reference_file,
        "output": output,
        "append": append,
        "omopolymeric_file": omopolymeric_file,
        "create_omopolymeric_file": create_omopolymeric_file,
        "splicing_file": splicing_file,
        "remove_header": args.remove_header,
        "bed_file": bed_file
        }
    
#     print("RUNNING REDItools 2.0 with the following options", options)
    
    return options

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

    options = parse_options()
    
    analyze(options)
    
    
    
