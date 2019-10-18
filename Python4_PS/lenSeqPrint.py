#!/usr/bin/env python3
import sys

seq = ['ATGCCCGGCCCGGC','GCGTGCTAGCAATACGATAAACCGG', 'ATATATATCGAT','ATGGGCCC']
for s in seq:
	print(seq.index(s), '\t', len(s), '\t', s)
print('Done!')

temp = []
[temp.append((len(s), s)) for s in seq]
print(temp)
