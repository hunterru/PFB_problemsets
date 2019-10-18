#!/usr/bin/env Python3
import math

nums = [101,2,15,22,95,33,2,27,72,15,52]

for num in nums :
	if num%2 == 0 :
		print(num)
print('Done!')

ev = 0
od = 0
for num in nums :
	print(num)
	if num%2 == 0 :
		ev += num
	else :
		od += num

print('Sum of even number:' , ev)
print('Sum of adds:' , od)
		
