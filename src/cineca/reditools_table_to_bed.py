import sys
import os
import gzip

import argparse	
if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='REDItools 2.0 table to BED file converter')
	parser.add_argument('-i', '--table-file', help='The RNA-editing events table to be converted')
	parser.add_argument('-o', '--bed_file', help='The output bed file')
	args = parser.parse_known_args()[0]

	input = args.table_file
	output = args.bed_file
	
	input_root, ext = os.path.splitext(input)
	fd_input = gzip.open(input, "r") if ext == ".gz" else open(input, "r")
	fd_output = open(output, "w")
	
	LOG_INTERVAL = 1000000
	last_chr = None
	start = None
	end = None
	
	total = 0
	with fd_input:
		for line in fd_input:
			total += 1
			if total % LOG_INTERVAL == 0:
				sys.stderr.write("[{}] {} lines read from {}\n".format(last_chr, total, input))
				
			fields = line.strip().split()
			chr = fields[0]
			pos = int(fields[1])
			
			if last_chr != chr or (end is not None and pos > end + 1):
				if last_chr is not None and start is not None and end is not None:
					fd_output.write("{}\t{}\t{}\n".format(last_chr, start, end))
					start = None
					end = None
			
			if start is None:
				start = pos
			
			if last_chr != chr:
				last_chr = chr
			
			if end is None or pos == end + 1:
				end = pos
			
	if last_chr is not None and start is not None and end is not None:
		fd_output.write("{}\t{}\t{}\n".format(last_chr, start, end))
		start = None
		end = None
		
	sys.stderr.write("{} lines read from {}\n".format(total, input))
	