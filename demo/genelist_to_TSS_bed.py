#!/usr/bin/env python3


#################################################
# FUNCTION: Check if a line is a single string
#################################################
def singlestring(line):
	import re
	if len(re.findall(r"(\w+)",line)) > 1:
		return False
	else:
		return True

#################################################################
# FUNCTION: Make a dictionary from a bed file (keys are ENS IDs)
#################################################################
def bedfile_to_dict_ENS(bedfile):
	gene_dict = {}

	# Open bed file, save all the values to a dictionary
	with open(bedfile, 'r') as fo_bed:

		for line in fo_bed:
			line = line.rstrip()
			
			field_str = "chr start end namefield score strand"
			fields = field_str.split(' ')
			genedata = dict(zip(fields,line.split('\t')))
			
			# Get the Ensembl ID name, use as IDs for dictionary
			namefield = genedata.pop('namefield')
			namefield_regex = re.search(r"^ID=gene:(\w+);Name=(\w+)" , namefield)	
			name = namefield_regex.group(1)
			genedata["name_other"] = namefield_regex.group(2)
			gene_dict[name] = genedata

		return gene_dict

##############################################################################
# FUNCTION: Make a dictionary from a bed file (keys are gene names)
##############################################################################
def bedfile_to_dict_genename(bedfile):
	
	import re
	gene_dict = {}

	with open(bedfile, 'r') as fo_bed: # open bed file and parse through

		for line in fo_bed:
			# Parse through line by line
			line = line.rstrip()
			field_str = "chr start end namefield score strand"
			fields = field_str.split(' ')
			genedata = dict(zip(fields,line.split('\t'))) # make dictionary of all info in bed line
			
			# Get the gene name type that you want, use as IDs for dictionary
			namefield = genedata.pop('namefield')
			namefield_regex = re.search(r"^ID=gene:(\w+);Name=(\w+)" , namefield)	
			name = namefield_regex.group(2)
			ENS_ID = namefield_regex.group(1)
			if name not in gene_dict:
				gene_dict[name] = {}
			gene_dict[name][ENS_ID] = genedata # each ENSEMBL ID is a separate subdictionary key
		
	return gene_dict

######################################################################################
# FUNCTION: Given a gene name and a gene dictionary, print out a bed file format line#
######################################################################################
def gene_to_bedline(gene, gene_dict):
	
	# if the gene is in the dictionary, return a bed-type line as a string
	if gene in gene_dict:
		# if the dictionary keys are ENS IDS
		if 'chr' in gene_dict[gene].keys():
			geneinfo = gene_dict[gene]
			bedline_to_print = "{}\t{}\t{}\t{} {}\t{}\t{}".format(geneinfo['chr'],geneinfo['start'],geneinfo['start'],gene,geneinfo['name_other'],geneinfo['score'],geneinfo['strand'])
			return bedline_to_print
		# if the dictionary keys are gene names
		else:
			bedline_to_print = ""
			ENS_ID_list = [] # list of ENS IDs mapping to the gene
			for ENS_ID in gene_dict[gene].keys():
				geneinfo = gene_dict[gene][ENS_ID]
				bedline_to_print += "{}\t{}\t{}\t{} {}\t{}\t{}\n".format(geneinfo['chr'],geneinfo['start'],geneinfo['start'],ENS_ID,gene,geneinfo['score'],geneinfo['strand'])
				ENS_ID_list.append(ENS_ID)
			# Provide note to user if the genename maps to multiple Ensembl IDs
			if len(ENS_ID_list) > 1:
				print("\nNote: The genename {} maps to multiple ENSIDs: {}\n".format(gene, " ".join(ENS_ID_list)))

			return bedline_to_print.rstrip()
	# if the gene is not in the dictionary, return false
	else:
		return False

############################################################################################
# FUNCTION: Given a TSS bed file and a list of genes, output only TSSs of genes of interest
############################################################################################
def TSSs_of_interest_bed(TSSbedfile, genelistfile, outputfolder="."):

	import os

	# Make a list of genes from all the user-provided info
	genelist = []
	# Read genes in the file if the user provided a file
	if genelistfile.endswith(".txt"):
		with open(genelistfile, 'r') as fo_genelist:
			for rawline in fo_genelist:
				line = rawline.rstrip()
				if not singlestring(line): # Quit if line of file contains more than one string
					print('\n!!!\nERROR! Your genelist file is not formatted correctly. Each line should contain one gene name!\n!!!\n')
					exit(1)
				genelist.append(line)
	# If it's not a file, then read in the genelistfile
	else:
		genelist.append(genelistfile)

	# Check if the genelist has Ensembl gene IDs to determine what "type" of genelist you have
	# Make the appropriate type of gene dictionary
	if genelist[0].startswith("ENS"):
		gene_dict = bedfile_to_dict_ENS(TSSbedfile)
	else:
		gene_dict = bedfile_to_dict_genename(TSSbedfile)

	# Go through the list of genes, make a new bed file with the TSSs of interest
	outputfilename = outputfolder+"/"+os.path.splitext(os.path.basename(TSSbedfile))[0]+".genesofinterest.bed"

	unmappedlist = [] # list of unmapped genes
	with open(outputfilename, 'w') as fo_output:
		for gene in genelist:
			bedline_to_print = gene_to_bedline(gene, gene_dict)
			if bedline_to_print:
				fo_output.write(bedline_to_print+"\n")
			else:
				unmappedlist.append(gene)
	print("\nThe TSSs that match your genes of interest were printed to the following bedfile:\n\t{}\n".format(outputfilename)) 
	
	# Print note to user if there were any unmapped gene names / IDs
	# Print out these genes to a new file so user can check the unmapped names
	unmappedcount = len(unmappedlist)
	if unmappedcount > 0:
		unmappedgenesfile = outputfolder+"/"+os.path.splitext(os.path.basename(genelistfile))[0]+".unmappedgenes.txt"
		with open(unmappedgenesfile, 'w') as fo_unmapped:
			for unmappedgene in unmappedlist:
				fo_unmapped.write(gene+"\n")
		print("\nNOTE: Your genelist had {} unmapped genes.\nThese genes are listed in {}.\n".format(unmappedcount, unmappedgenesfile))

	# Return outputfilename
	return outputfilename

#####################
# Main function
#####################

def main():
	
	from gff3_to_TSSbed import gff3_to_TSSbed
	import sys
	import os

	if len(sys.argv) < 3:
		usagemessage = "\t{:20} {} annotationfile.gff3 {}"
		print("\n\t!!!PROGRAM USAGE ERROR!!!\n")
		print(usagemessage.format("For a single gene:",sys.argv[0], "genename"))
		print(usagemessage.format("For a gene list:",sys.argv[0], "genelist.txt"))
		print("\n")
		sys.exit(1)

	gff3file = sys.argv[1]
	genelistfile = sys.argv[2]

	gff3_to_TSSbed(gff3file)

	TSSbedfile = os.path.splitext(gff3file)[0]+".TSS.bed"

	# Call function to print out a bed file from a bed file with a TSS of interest
	# and with a list of genes of interest
	TSSs_of_interest_bed(TSSbedfile, genelistfile) # can use a .txt file with a gene list

if __name__ == "__main__":
	main()
