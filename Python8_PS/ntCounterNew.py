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

seqs = {}
for geneName in fastaDic:
	geneName = geneName.rstrip()
	if geneName not in seqs:
		seqs[geneName] = {}
	if 'A' not in seqs[geneName]:
			seqs[geneName]['A'] = fastaDic[geneName].count('A')
	if 'C' not in seqs[geneName]:
			seqs[geneName]['C'] = fastaDic[geneName].count('C')
	if 'T' not in seqs[geneName]:
			seqs[geneName]['T'] = fastaDic[geneName].count('T')
	if 'G' not in seqs[geneName]:
			seqs[geneName]['G'] = fastaDic[geneName].count('G')


for Name in seqs:
	for base in seqs[Name]:
		out = seqs[Name][base]
		print(Name,base,out,sep='\t')
