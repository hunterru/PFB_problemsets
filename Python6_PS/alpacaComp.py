#!/usr/bin/env python3

allGeneList = []
with open('alpaca_all_genes.tsv','r') as allGenes:
	for lines in allGenes:
		temp = lines.split()
		allGeneList.append(temp[0])	
allGeneList = allGeneList[1:]
allGeneList = set(allGeneList)

allStemList = []
with open('alpaca_stemcellproliferation_genes.tsv','r') as stemGenes:
	for lines in stemGenes:
		temp = lines.split()
		allStemList.append(temp[0])
allStemList = allStemList[1:]
allStemList = set(allStemList)
	
allPigList = []
with open('alpaca_pigmentation_genes.tsv','r') as pigGenes:
	for lines in pigGenes:
		temp = lines.split()
		allPigList.append(temp[0])
allPigList = allPigList[1:]
allPigList = set(allPigList)

noPro = allGeneList - allStemList
print('The number of total genes is:', len(allGeneList))
print('The number of proliferation genes is:', len(allStemList))
print('The number of pigmentation genes is:', len(allPigList))
print('The number of genes not related to proliferation is:', len(noPro))

proAndPig = allStemList & allPigList
print('The number of genes related to BOTH proliferation and pigmentation is:', len(proAndPig))
