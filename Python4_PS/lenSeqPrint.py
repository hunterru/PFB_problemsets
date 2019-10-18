#!/usr/bin/env python3
import sys

seq = ['ATGCCCGGCCCGGC','GCGTGCTAGCAATACGATAAACCGG', 'ATATATATCGAT','ATGGGCCC']
count = 0
for s in seq:
	print(count, '\t', len(s), '\t', s)
	count += 1
print('Done!')

temp = []
[temp.append((len(s), s)) for s in seq]
print(temp)
