#!/usr/bin/env python3
import re

ID = []
with open('Python_08.fasta') as sequence:
	for line in sequence:
		line = line.rstrip()
		out = re.search('^>\S+' , line)
		if out:
			print(type(out))
