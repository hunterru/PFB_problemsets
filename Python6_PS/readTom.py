#!/usr/bin/env python3

with open('Python_06.txt','r') as song, open('Python_06_uc.txt','w') as songWrite:
	for line in song:
		line = line.rstrip()
		line_up = line.upper()
		songWrite.write(line_up+'\n')
		print(line_up)
