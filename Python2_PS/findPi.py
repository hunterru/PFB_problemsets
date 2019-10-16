#!/usr/bin/env python3
import sys	# to read input from command line
import math	# to use advanced math functions

sampleNum = sys.argv[1]
sampleNum = float(sampleNum)
PI = 3.14159

if sampleNum < PI :
	output1 = 'Unfortunately,'
	output2 = 'is less than pi :('
elif sampleNum == PI :
	output1 = 'You did it!'
	output2 = 'is pi! You found pi! :)'
else :
	output1 = 'Unfortunately,'
	output2 = 'is greater than pi :('
print(output1, sampleNum, output2)
