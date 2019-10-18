#!/usr/bin/env python3

countLine = 0
countChar = 0
with open('Python_06.fastq') as dataFastq:
	for line in dataFastq:
		countLine += 1
		countChar += len(line)
print('This is the total number of lines:',countLine)
print('This is the total number of characters:',countChar)
print('This is the average length of each line:',countChar/countLine)
