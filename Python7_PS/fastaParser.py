#!/usr/bin/env python3
import re

fastaDic = {}
with open('Python_07.seq.reverse.txt') as fastaData :
	for line in fastaData :
		out = re.search('^>(\S+)', line)
		if out :
			ID = out.group(1)
		else:
			fastaDic[ID] = line.rstrip()
print(fastaDic)
