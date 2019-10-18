#!/usr/bin/env python3
import sys

things = {
	'book' : 'Wheel of Time' ,
	'song' : 'Tom Petty - I Won\'t Back Down' ,
	'tree' : 'Aspen'
	}

things['organism'] = 'Humans'

# Using a string conditional to form a list of the keys
opt = [key for key in things.keys()]
opt_join = ', '.join(opt)

print('What value from the following key would you like to see? Key:' , opt_join)
userIn = input()

fav_thing = userIn
print(things[fav_thing])

# Changing the organism based on user input
print('Now, let\'s change the organism. What is your favorite organism?')
orgIn = input()
things['organism'] = orgIn
print('The following is your new dictionary')
for key in things:
	thing = things[key]
	print("{:>15}".format(key), "{:<5}".format(thing))
