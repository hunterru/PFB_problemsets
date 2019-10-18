#!/usr/bin/env Python3

taxa = 'sapiens, erectus, neaderthalensis'
print(taxa)
print(taxa[1])

species  = taxa.split(', ')
print(species)
print(species[1])
print(type(species))  # species is a 'list'

print(sorted(species))

print(sorted(species,key=len))
