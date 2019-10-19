#!/usr/bin/env python3
import re

newSeq = ''
with open('Python_07_ApoI.fasta') as sequence:
	for line in sequence:
		out = re.search('^>\S+' , line)
		if not out :
			temp = line.rstrip()
			newSeq = newSeq+temp

	info = re.findall('[AG]AATT[TC]' , newSeq)
	for seq in info:
		seqHat = seq[:1] + '^' +seq[1:]
		newSeq = newSeq.replace(seq,seqHat)
	print(newSeq)
	
	newSeqList = newSeq.split('^')
	print(type(newSeqList))
	print(sorted(newSeqList,key=len))
