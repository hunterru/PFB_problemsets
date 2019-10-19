#!/usr/bin/env python3
import re

count = 1
with open('Python_07_nobody.txt') as poem:
	for line in poem:
		for found in re.finditer(r"Nobody" , line):
			print(type(found))
			out = found.start(0)
			print('The line is:' , count , 'and the position is:' , out)
		count += 1
print('Done!')


with open('Python_07_nobody.txt','r') as poem, open('FavName.txt','w') as out:
   for line in poem:
      for found in re.sub(r"Nobody" , 'FavName' ,  line):
         out.write(found)
print('Done!')
