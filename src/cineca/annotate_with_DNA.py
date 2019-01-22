import sys
import os
import gzip

columns = {
	"Region": 0,
	"Position": 1,
	"Reference": 2,
	"Strand": 3,
	"Coverage": 4,
	"MeanQ": 5,
	"BaseCount": 6,
	"AllSubs": 7,
	"Frequency": 8,
	"gCoverage": 9,
	"gMeanQ": 10,
	"gBaseCount": 11,
	"gAllSubs": 12,
	"gFrequency": 13
}

def read_line(fd):
	line = next(fd, None)
	return line.strip().split("\t") if line is not None else []

def is_smaller_than_or_equal_to(fields1, fields2):
# 	if not fields1 and not fields2: return True
# 	if fields1 and not fields2: return True
	if not fields1 and fields2: return False
	
	if fields1 and fields2:
		region1 = get(fields1, "Region")
		region2 = get(fields2, "Region")
		
		index1 = chromosomes.index(region1) if region1 in chromosomes else chromosomes.index("chr" + region1)
		index2 = chromosomes.index(region2) if region2 in chromosomes else chromosomes.index("chr" + region2)
		
# 		sys.stderr.write(" ".join([str(x) for x in [region1, region2, index1, index2]]) + "\n")
		
		if index1 < index2:
			return True
		
		if index1 > index2:
			return False
		
		return index1 == index2 and int(get(fields1, "Position")) <= int(get(fields2, "Position"))
	
	return True

def get(fields, column):
	value = None
	
	index = columns[column]
	if len(fields) >= index:
		value = fields[index]
		 
	return value

comp = {'A':'T','T':'A','C':'G','G':'C'}
indexes = {v: k for k, v in dict(enumerate('ACGT')).iteritems()}

def annotate(fields1, fields2):
	
	strand = get(fields1, "Strand")
	
	if strand == '0':
		base_count = eval(get(fields2, "BaseCount")) # BaseCount[A,C,G,T]
		fields2[columns["BaseCount"]] = str([base_count[indexes[comp[b]]] for b in 'ACGT'])
		
		subs = get(fields2, "AllSubs").split(" ")
		fields2[columns["AllSubs"]] = " ".join([''.join([comp[b] if b != "-" else b for b in sub]) for sub in subs])
	
	for field in ["Coverage", "MeanQ", "BaseCount", "AllSubs",  "Frequency"]:
		annotation = get(fields2, field)
# 		if annotation is None:
# 			print(fields1)
# 			print(fields2)
# 			print(field, annotation)
		
		fields1[columns["g" + field]] = annotation

chromosomes = []
def load_chromosomes(fai):
	with open(fai, "r") as reader:
		for line in reader:
			chromosome = line.strip().split("\t")[0]
			if chromosome in chromosomes: continue
			chromosomes.append(chromosome)

LOG_INTERVAL = 1000000

import argparse	
if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='REDItools 2.0 annotator')
	parser.add_argument('-r', '--rna-file', required=True, help='The RNA-editing events table to be annotated')
	parser.add_argument('-d', '--dna-file', required=True, help='The RNA-editing events table as obtained from DNA-Seq data')
	parser.add_argument('-R', '--reference', required=True, help='The .fai file of the reference genome containing the ordered chromosomes')
	parser.add_argument('-Z', '--only-omozygotes', default=False, action='store_true', help='Exclude positions with multiple changes in DNA-Seq')
	args = parser.parse_known_args()[0]

	file1 = args.rna_file
	file2 = args.dna_file
	fai_file = args.reference
	load_chromosomes(fai_file)
	only_omozygotes = args.only_omozygotes
	
	sys.stderr.write("[INFO] {} CHROMOSOMES LOADED\n".format(len(chromosomes)))
	
	file1root, ext1 = os.path.splitext(file1)
	file2root, ext2 = os.path.splitext(file2)
	
	fd1 = gzip.open(file1, "r") if ext1 == ".gz" else open(file1, "r")
	fd2 = gzip.open(file2, "r") if ext2 == ".gz" else open(file2, "r")
	fd3 = sys.stdout
	
	total1 = 0
	total2 = 0
	last_chr = None
	with fd1, fd2, fd3:
	
		fields1 = read_line(fd1)
		total1 += 1
		if fields1[0] == "Region":
			fields1 = read_line(fd1)
			total1 += 1
		
		fields2 = read_line(fd2)
		total2 += 1
		if fields2[0] == "Region":
			fields2 = read_line(fd2)
			total2 += 1

		while fields1 or fields2:
			
			if fields1[0] != last_chr:
				last_chr = fields1[0]
				sys.stderr.write("ANALYZING CHROMOSOME " + last_chr + "\n")
			
			f1_less_than_f2 = is_smaller_than_or_equal_to(fields1, fields2)
			f2_less_than_f1 = is_smaller_than_or_equal_to(fields2, fields1) 
			are_equal = f1_less_than_f2 and f2_less_than_f1

# 			sys.stderr.write(str(fields1) + "\n")
# 			sys.stderr.write(str(fields2) + "\n")
# 			sys.stderr.write(str(f1_less_than_f2) + " " + str(f2_less_than_f1) + " " + str(are_equal) + "\n")
# 			raw_input()

			omozigote = True if not fields2 else not are_equal or fields2[columns["AllSubs"]] == "-"
			
			if are_equal:
				annotate(fields1, fields2)
			
			if fields1:
				if not only_omozygotes or omozigote:
					fd3.write("\t".join(fields1) + "\n")
				else:
					sys.stderr.write("[INFO] [{}] Discarding {}:{} because DNA data is not omozygote from {}\n".format(last_chr, fields1[0], fields1[1], file1))
				
			if f1_less_than_f2:
				fields1 = read_line(fd1)
				total1 += 1
				
			if f2_less_than_f1:
				fields2 = read_line(fd2)
				total2 += 1
				
			if total1 % LOG_INTERVAL == 0:
				sys.stderr.write("[INFO] [{}] {} lines read from {}\n".format(last_chr, total1, file1))
				
			if total2 % LOG_INTERVAL == 0:
				sys.stderr.write("[INFO] [{}] {} lines read from {}\n".format(last_chr, total2, file2))
		
	sys.stderr.write("[INFO] {} lines read from {}\n".format(total1, file1))
	sys.stderr.write("[INFO] {} lines read from {}\n".format(total2, file2))
	