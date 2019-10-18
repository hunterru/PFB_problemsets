#!/usr/bin/env python3

allStemList = []
with open('alpaca_stemcellproliferation_genes.tsv','r') as stemGenes:
   for lines in stemGenes:
      temp = lines.split()
      allStemList.append(temp[0])
allStemList = allStemList[1:]
allStemList = set(allStemList)

allTransList = []
with open('alpaca_transcriptionFactors.tsv','r') as transGenes:
   for lines in transGenes:
      temp = lines.split()
      allTransList.append(temp[0])
allTransList = allTransList[1:]
allTransList = set(allTransList)

transPro = allStemList & allTransList
print('The number of proliferation genes is:', len(allStemList))
print('The number of transcription genes is:', len(allTransList))
print('The number of transcription factors for proliferation is:', len(transPro))
