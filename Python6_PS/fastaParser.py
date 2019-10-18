#!/usr/bin/env python3

fastaDic = {}
with open('Python_06.seq.reverse.txt') as fastaData :
	for line in fastaData :
		if line[0] == '>':
			temp = line[1:]
			name,junk = temp.split('\t')
		else:
			fastaDic[name] = line
print(fastaDic)
