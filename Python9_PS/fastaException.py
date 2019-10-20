#!/usr/bin/env python3
import re
import sys

file = ''
try:
	file = sys.argv[1]
	print('User provided file name:' , file)
	if not file.endswith('.fa') or file.endswith('.fasta') or file.endswith('.nt'):
		raise ValueError('Not a FASTA file (idiot!)')
	fastaDic = {}
	with open(file) as fastaData :
		for line in fastaData :
			out = re.search('^>(\S+)', line)
			if out :
				ID = out.group(1)
			else:
				fastaDic[ID] = line.rstrip()
				temp = line.strip()
				temp = temp.upper()
				for nt in temp:
					nucs = 'ATCGN'
					if not nucs.count(nt): 
						print(nt)
						raise ValueError('An unaccepted nucleotide value is present in the sequence:', nt)
				
	print('Try complete!')
except IndexError:
	print('Please provide a file name')
except IOError:
	print('Can\'t find file:' , file)

