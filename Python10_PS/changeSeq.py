#!/usr/bin/env python3

def dna_lengthCut(ins,clipLength):
	DNA = ins
	cutter = clipLength
	print('This is the original sequence:\n'+DNA+'\n')
	dna =''
	for line in DNA:
		line = line.rstrip()
		dna = dna+line

	count = 0	
	newSeq = ''
	for nt in dna:
		count += 1
		if count%cutter != 0:
			newSeq = newSeq+nt
		else:
			newSeq = newSeq+nt+'\n'
	print('This is the modified sequence:\n'+newSeq+'\n')


def gc_content(ins):
	DNA = ins
	nucs = 'ACTG'
	gCount = DNA.count('G')
	cCount = DNA.count('C')
	temp = (gCount + cCount) / len(DNA)
	return('{:.2%}'.format(temp))

def revComp(ins):
	DNA = ins
	DNA = DNA.upper()
	DNA = DNA.replace('A','a')
	DNA = DNA.replace('C','c')
	DNA = DNA.replace('a','t')
	DNA = DNA.replace('c','g')
	DNA = DNA.replace('G','C')
	DNA = DNA.replace('T','A')
	DNA = DNA.upper()
	return(DNA)



DNA1 = 'GATGGGATTGGGGTTTTCCCCTCCCATGTGCTCAAGACTGGCGCTAAAAGTTTTGAGCTTCTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCT'

DNA2 = '''GATGGGATTGGGGTTTTCCCCTCCCATGTGCTCAAGACTGGCGCTAAAAGTTTTGAGCTTCTCAAAAGTCTAGAGCCACC
GTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACG
CTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCC
TCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAA
TGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATG
CCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCT
GTCATCTTCT'''

#dna_lengthCut(DNA2,70)
#gc_out = gc_content(DNA1)
revOut = revComp(DNA1)
print(revOut)