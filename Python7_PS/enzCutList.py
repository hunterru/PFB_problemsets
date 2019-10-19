#!/usr/bin/env python3
import re

with open('bionet.txt','r') as textFile:
	text = textFile.readlines()[10:]
	for line in text:
		temp = re.search(r'(\S+\s{1,3}\S+)', line)
		enz = temp.group(1)
		print(enz)
		
