#!/usr/bin/env python3
import re

with open('Python_08.fasta') as sequence:
	fastaDic = {}
	for line in sequence:
		line = line.rstrip()
		out = re.search('^>(\S+)' , line)
		if out:
			ID = out.group(1)
			seq = ''
		else:
			seq += line
			fastaDic[ID] = seq

with open('Python_08.codons-6frame.nt','w') as codonOutput:
	for geneName in fastaDic:
		geneName = geneName.rstrip()
		codonOutput.write(geneName+'-frame-1-codons\n')
		seq = fastaDic[geneName]
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '
				count = 0
		codonOutput.write(seqSpace+'\n')
		
		codonOutput.write(geneName+'-frame-2-codons\n')
		seq = fastaDic[geneName]
		seq = seq[1:]
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '
				count = 0
		codonOutput.write(seqSpace+'\n')
		
		codonOutput.write(geneName+'-frame-3-codons\n')
		seq = fastaDic[geneName]
		seq = seq[2:]
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '
				count = 0
		codonOutput.write(seqSpace+'\n')

*** you are here - need to reverse complement!	
	for geneName in fastaDic:
		geneName = geneName.rstrip()
		codonOutput.write(geneName+'-frame-1-codons\n')
		seq = fastaDic[geneName]
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '
				count = 0
		codonOutput.write(seqSpace+'\n')
		
		codonOutput.write(geneName+'-frame-2-codons\n')
		seq = fastaDic[geneName]
		seq = seq[1:]
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '
				count = 0
		codonOutput.write(seqSpace+'\n')
		
		codonOutput.write(geneName+'-frame-3-codons\n')
		seq = fastaDic[geneName]
		seq = seq[2:]
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '
				count = 0
		codonOutput.write(seqSpace+'\n')
