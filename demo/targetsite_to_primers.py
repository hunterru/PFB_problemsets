#!/usr/bin/env python3

import os
import re
import sys

#########################################
# FUNCTION: Return revcomp of a sequence
#########################################
def revcomp(seq):
	seq = seq.upper()
	seq_comp = seq.replace("A","t")
	seq_comp = seq_comp.replace("T","a")
	seq_comp = seq_comp.replace("C","g")
	seq_comp = seq_comp.replace("G","c")
	seq_revcomp = seq_comp[::-1].upper()
	return seq_revcomp

################################################
# FUNCTION: Add G to first pos if there is no G
################################################
def startG(seq):
	if seq.startswith("G"):
		return(seq)
	else:
		return("G"+seq)

###############################################################
# FUNCTION: Returns primer parameters for a given strategy
###############################################################
def cloning_parameters(cloningstrategy):
	
	import os

	# Find the primer parameters for given cloning strategy from a file
	file_with_cloning_info = "cloningstrategy_primers.txt"
	if not os.path.isfile(file_with_cloning_info):
		print("\n\tMissing file: 'cloningstrategy_primers.txt'\n")

	with open(file_with_cloning_info,'r') as fo:
		for line in fo:
			line = line.rstrip()
			strategy_info = line.split("\t")
			if cloningstrategy == strategy_info[0]:
				Grequired = strategy_info[1]
				p1specs = strategy_info[2]
				p2specs = strategy_info[3]
				break

	# Return the parameters. If couldn't find cloning strategy in file, provide helpful error message
	try:
		return (Grequired, p1specs, p2specs)
	except NameError:
		print("\n\tError! Could not find cloning strategy {} in the file 'cloningstrategy_primers.txt.'\n".format(cloningstrategy))

########################################################
# FUNCTION: Returns primers to order given a target site
#########################################################
def targetsite_to_primers(seq, cloningstrategy, Grequired, p1specs, p2specs):

	# Add G to seq if a G is required
	if Grequired == "G":
		seq = startG(seq)
	
	# Get revcomp of the seq
	seq_revcomp = revcomp(seq)

	# Get the sequence of primer1
	p1specs_list = p1specs.lower().split(" ")
	if len(p1specs_list) == 2:
		primer1 = p1specs_list[0] + seq + p1specs_list[1]
	else:
		primer1 = p1specs_list[0] + seq

	# Get the sequence of primer2
	if "IVT" in cloningstrategy:
		primer2 = p2specs.lower()
	else:
		p2specs_list = p2specs.lower().split(" ")
		if len(p2specs_list) == 2:
			primer2 = p2specs_list[0] + seq_revcomp + p2specs_list[1]
		else:
			primer2 = p2specs_list[0] + seq_revcomp

	# Return primer sequences
	return (primer1, primer2)

##########################################
# MAIN                                   
##########################################
def main():
	
	# Check for 2 arguments. Raise error message and exit if not enough arguments.
	import sys
	if len(sys.argv) < 3:
		print("\n\tPROGRAM USAGE:\n")
		print("\t{} sequence cloningstrategy".format(sys.argv[0]))
		print("\tThe cloning strategy should be listed in 'cloningstrategy_primers.txt'\n")
		sys.exit(1)

	# Get arguments from the command line
	seq = sys.argv[1]
	cloningstrategy = sys.argv[2]
	
	# Run commands	
	(Grequired, p1specs, p2specs) = cloning_parameters(cloningstrategy)

	(primer1, primer2) = targetsite_to_primers(seq, cloningstrategy, Grequired, p1specs, p2specs)
	print("Primer1:",primer1)
	print("Primer2:",primer2)



if __name__ == "__main__":
	main()
