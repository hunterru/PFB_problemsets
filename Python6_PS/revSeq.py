#!/usr/bin/env python3

with open('Python_06.seq.txt','r') as sequence:
	for line in sequence:	
		marker,seq = line.split()
		revLine = seq[::-1]
		newRevLine = ''
		for pos in seq:
			if pos == 'A':
				newPos = 'T'
			elif pos == 'T':
				newPos = 'A'
			elif pos == 'C':
				newPos = 'G'
			else:
				newPos = 'C'
			newRevLine = newRevLine+newPos
		print('>'+marker+'\tthis is the reverse complement\n',end='')
		print(newRevLine+'\n',end='')
