#!/usr/bin/env python3
import re

enzDict = dict()
with open('bionet.txt','r') as textFile:
	text = textFile.readlines()[10:]
	for line in text:
		found = line.split(' ')
		if found:
			enz = found[0]+' '+found[1]
			enz = enz.rstrip()
			site = found[-1]
			site = site.rstrip()
			enzDict[enz] = site
print(enzDict)
