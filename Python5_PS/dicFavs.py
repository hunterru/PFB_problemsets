#!/usr/bin/env python3
import sys

things = {
	'book' : 'Wheel of Time' ,
	'song' : 'Tom Petty - I Won\'t Back Down' ,
	'tree' : 'Aspen'
	}

things['organism'] = 'Humans'

print([key for key in things.keys()])

print('What value from the following key would you like to see? Key:' , things.keys() )
userIn = input()

fav_thing = userIn
print(things[fav_thing])
