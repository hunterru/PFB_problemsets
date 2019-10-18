#!/usr/bin/env python3
import sys  # to read input from command line

begin = int(sys.argv[1])
end = int(sys.argv[2])
end = end + 1
type(begin)
type(end)

rg  = range(begin, end)
for num in rg :
	if num%2 != 0 :
		print(num)
print('Done!')
