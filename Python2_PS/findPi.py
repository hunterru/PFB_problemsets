#!/usr/bin/env python3
import sys	# to read input from command line
import math	# to use advanced math functions

print('Would you like to play a game (Yes/No)?')
ans = input()
if ans == 'Yes' :
	print('What is your best guess for pi?')	
	sampleNum  = input()	
		#	sampleNum = sys.argv[1]
	sampleNum = float(sampleNum)
	pi = 3.14
	precPoor = 0.1
	precGood = 0.01

	if round(abs(pi - sampleNum),2) > precPoor :
		output1 = 'Not even close!'
		output2 = 'is really far from pi!'
		s = 'Your precision is off by > |{n}|!'
		output3 = s.format(n=precPoor)
	elif round(abs(pi - sampleNum),2) >=  precGood :
		output1 = 'Not too bad.'
		output2 = 'is not quite equal to pi.'
		s = 'Your precision is off by >= |{n}|.'
		output3 = s.format(n=precGood)
	elif round(abs(pi - sampleNum),2) < precGood :
		output1 = 'You did it!'
		output2 = 'is pi! You found pi! :)'
		s = 'Your precision is < |{n}|.'
		output3 = s.format(n=precGood)
	print(output1, sampleNum, output2, output3)
else :
	print('You are lame! Who doesn\'t like a game!')
