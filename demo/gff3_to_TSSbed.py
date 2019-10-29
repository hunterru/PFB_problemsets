#!/usr/bin/env python3

import sys
import os
import re

# Function that takes gff3filename as input
# Outputs a bed file for the TSS

def gff3_to_TSSbed(gff3filename, outputfolder="."):

	
	# Determine the output file name
	outputfilename = outputfolder+"/"+os.path.splitext(gff3filename)[0]+".TSS.bed"

	# fields in a gff3 files
	field_str = "seqid source genetype start end score strand phase attributes"
	fields = field_str.split(' ')
	
	# open gff3 file, read in lines with gene information
	with open(gff3filename, 'r') as fo, open(outputfilename, 'w') as outputfile:
		
		for line in fo:
			line = line.strip()
			
			# skip header lines
			if line.startswith('#'):
				continue
			# get info from non-header lines
			else:
				data = dict(zip(fields,line.split('\t')))		
				# only keep lines that are of type "gene"
				# print relevant columns to file

				if data['genetype'] == 'gene':

					attributes = data['attributes'].split(";")
					gene_name = attributes[0]+";"+attributes[1]
					strand = data['strand']
										
					bed_line = "chr{}\t{}\t{}\t{}\t{}\t{}\n"
					outputfile.write(bed_line.format(data['seqid'], data['start'], data['start'], gene_name, '', strand))
	
	# Tell user where the TSS file was outputted
	print("\nTSSs matching all genes were printed to the following file in bedfile format:\n\t{}\n".format(outputfilename))
	return outputfilename

def main():

	import sys
	if len(sys.argv) < 2:
		print("\n\t!!!PROGRAM USAGE ERROR!!!\n")
		print("\t{} annotationfile.gff3\n".format(sys.argv[0]))
		sys.exit(1)
	
	gff3filename = sys.argv[1]
	gff3_to_TSSbed(gff3filename)
	
if __name__ == "__main__":
	main()
