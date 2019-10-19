#!/usr/bin/env python3
import re

with open('Python_07.fasta') as data:
	for line in data:
		out = re.search(r"^(>\S+)(\s+.+)", line)
		if out:
			ident = out.group(1)
			comment = out.group(2)
			print('ID:' , ident , 'Description:' , comment)
print('Done!')
