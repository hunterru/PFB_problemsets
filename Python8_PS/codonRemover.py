#!/usr/bin/env python3
import re

translation_table = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}

with open('Python_08.fasta') as sequence:
	fastaDic = {}
	for line in sequence:
		line = line.rstrip()
		out = re.search('^>(\S+)' , line)
		if out:
			ID = out.group(1)
			seq = ''
		else:
			seq += line
			fastaDic[ID] = seq

with open('Python_08.codons-6frame.nt','w') as codonOutput, open('Python_08.translated.aa','w') as aminoAcids, open('Python_08.translated-longest.aa','w') as longAcid, open ('Python_08.orf-longest.nt','w') as longCodon:
	for geneName in fastaDic:
		longPeps = []
		longCods = []
		startPeps = []
		endPeps = []
		longSeqs = []
		geneName = geneName.rstrip()
		codonOutput.write(geneName+'-frame-1-codons\n')
		seq = fastaDic[geneName]
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '+nt
				count = 1
		codonOutput.write(seqSpace+'\n')
		
		aminoAcids.write(geneName+'-frame-1-proteins\n')
		codonParse = seqSpace.split()
		aa = ''
		for codon in codonParse:
			if len(codon) == 3:
				temp = translation_table[codon]
				aa += temp
		aminoAcids.write(aa+'\n')
		
		if re.findall(r'M\w+\*',aa):
			pepTemp = re.search(r'(M)\w+(\*)',aa)
			longPeps.append(pepTemp[0])
			startPeps.append(pepTemp.start(1))
			endPeps.append(pepTemp.start(2))
			longSeqs.append(seq)

## Starting Frame 2 - codon and AA		
		codonOutput.write(geneName+'-frame-2-codons\n')
		seq = fastaDic[geneName]
		seq = seq[1:]
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '+nt
				count = 1
		codonOutput.write(seqSpace+'\n')

		aminoAcids.write(geneName+'-frame-2-proteins\n')
		codonParse = seqSpace.split()
		aa = ''
		for codon in codonParse:
			if len(codon) == 3:
				temp = translation_table[codon]
				aa += temp
		aminoAcids.write(aa+'\n')

		if re.findall(r'M\w+\*',aa):
			pepTemp = re.search(r'(M)\w+(\*)',aa)
			longPeps.append(pepTemp[0])
			startPeps.append(pepTemp.start(1))
			endPeps.append(pepTemp.start(2))
			longSeqs.append(seq)

## Starting Frame 3 - codon and AA
		codonOutput.write(geneName+'-frame-3-codons\n')
		seq = fastaDic[geneName]
		seq = seq[2:]
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '+nt
				count = 1
		codonOutput.write(seqSpace+'\n')

		aminoAcids.write(geneName+'-frame-3-proteins\n')
		codonParse = seqSpace.split()
		aa = ''
		for codon in codonParse:
			if len(codon) == 3:
				temp = translation_table[codon]
				aa += temp
		aminoAcids.write(aa+'\n')


		if re.findall(r'M\w+\*',aa):
			pepTemp = re.search(r'(M)\w+(\*)',aa)
			longPeps.append(pepTemp[0])
			startPeps.append(pepTemp.start(1))
			endPeps.append(pepTemp.start(2))
			longSeqs.append(seq)

# now reverse complement for all 3 codon positions
## Starting Frame 4 - codon and AA
		codonOutput.write(geneName+'-frame-4-codons\n')
		seq = fastaDic[geneName]
		newSeq = seq[::-1]
		newSeq = newSeq.replace('A','a')
		newSeq = newSeq.replace('C','c')
		newSeq = newSeq.replace('a','t')
		newSeq = newSeq.replace('c','g')
		newSeq = newSeq.replace('T','A')
		newSeq = newSeq.replace('G','C')
		revSeq = newSeq.upper()
		
		seq = revSeq
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '+nt
				count = 1
		codonOutput.write(seqSpace+'\n')
		
		aminoAcids.write(geneName+'-frame-4-proteins\n')
		codonParse = seqSpace.split()
		aa = ''
		for codon in codonParse:
			if len(codon) == 3:
				temp = translation_table[codon]
				aa += temp
		aminoAcids.write(aa+'\n')


		if re.findall(r'M\w+\*',aa):
			pepTemp = re.search(r'(M)\w+(\*)',aa)
			longPeps.append(pepTemp[0])
			startPeps.append(pepTemp.start(1))
			endPeps.append(pepTemp.start(2))
			longSeqs.append(seq)

## Starting Frame 5 - codon and AA
		codonOutput.write(geneName+'-frame-5-codons\n')
		seq = revSeq[1:]
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '+nt
				count = 1
		codonOutput.write(seqSpace+'\n')
		
		aminoAcids.write(geneName+'-frame-5-proteins\n')
		codonParse = seqSpace.split()
		aa = ''
		for codon in codonParse:
			if len(codon) == 3:
				temp = translation_table[codon]
				aa += temp
		aminoAcids.write(aa+'\n')


		if re.findall(r'M\w+\*',aa):
			pepTemp = re.search(r'(M)\w+(\*)',aa)
			longPeps.append(pepTemp[0])
			startPeps.append(pepTemp.start(1))
			endPeps.append(pepTemp.start(2))
			longSeqs.append(seq)

## Starting Frame 6 - codon and AA
		codonOutput.write(geneName+'-frame-6-codons\n')
		seq = revSeq[2:]
		seqSpace = ''
		count = 0
		for nt in seq:
			if count < 3:
				seqSpace += nt
				count += 1
			else:
				seqSpace += ' '+nt
				count = 1
		codonOutput.write(seqSpace+'\n')

		aminoAcids.write(geneName+'-frame-6-proteins\n')
		codonParse = seqSpace.split()
		aa = ''
		for codon in codonParse:
			if len(codon) == 3:
				temp = translation_table[codon]
				aa += temp
		aminoAcids.write(aa+'\n')


		if re.findall(r'M\w+\*',aa):
			pepTemp = re.search(r'(M)\w+(\*)',aa)
			longPeps.append(pepTemp[0])
			startPeps.append(pepTemp.start(1))
			endPeps.append(pepTemp.start(2))
			longSeqs.append(seq)


		if longPeps:
			le = max(len(peps) for peps in longPeps)
			pepCountNew = 0
			for peps in longPeps:
				pepCountNew += 1
				if len(peps) == le:
					pepIndex = longPeps.index(peps)
					longAcid.write(geneName+'-frame-'+str(pepCountNew)+'-protein\n')
					longAcid.write(peps+'\n')
					longCodon.write(geneName+'-frame'+str(pepCountNew)+'-codon\n')
					startCodon = startPeps[pepIndex]*3
					endCodon = endPeps[pepIndex]*3 + 3
					codSeq = longSeqs[pepIndex]
					longCodon.write(codSeq[startCodon:endCodon]+'\n')
		else:
			longAcid.write(geneName+'-frame-Empty-protein\n')
			longAcid.write(' \n')
			longAcid.write(' \n')
			longCodon.write(geneName+'-frame-Empty-codon\n')
			longCodon.write(' \n')
			longCodon.write(' \n')		
