#!/usr/bin/env python3

import argparse
import os, sys, re, glob, shutil
import subprocess
import gzip
from gff3_to_TSSbed import gff3_to_TSSbed
from genelist_to_TSS_bed import TSSs_of_interest_bed
from overlap import range_overlap
from grna_gen import denovoGuideRnaAnno

##########-Start-###################
# Function: Download Gff3 file from Ensembl .gff3 ftp site

def downloadGff3(speciesName,outputDir):
	fileList = 'wget ' + '\'' + 'ftp://ftp.ensembl.org/pub/release-98/gff3/' + speciesName + '/CHECKSUMS' + '\''
	fileList_run = subprocess.run(fileList, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if fileList_run.returncode != 0:
		print("FAILED! Unable to locate list file. Check for existence on the Ensembl ftp site.")
		exit(2)
	else:
		with open('CHECKSUMS') as gZList:
			for lines in gZList:
				lines = lines.rstrip()
				if lines.endswith('.98.gff3.gz'):
					temp = lines.split()
					fileGrab = temp[2]

		siteCmd = 'wget ' + '\'' + 'ftp://ftp.ensembl.org/pub/release-98/gff3/' + speciesName + '/' + fileGrab + '\''
		siteCmd_run = subprocess.run(siteCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if siteCmd_run.returncode != 0:
			print("FAILED! Either your species is not available or the species did not match. Please check the Ensembl FTP site and try again")
			exit(2)
		else:
			contents=os.listdir(outputDir)
			for files in contents:
				if files.endswith('98.gff3.gz'):
					preGz = files           # the file name that matches the desired gff3 file
					with gzip.open(preGz,'rb') as f_in:
						with open(preGz[:-3], 'wb') as f_out:
							shutil.copyfileobj(f_in, f_out)  # unzipping the gff3 file
			newName = preGz[:-3]    # removing the '.gz' from the name since it is now unzipped
			os.remove(preGz)        # deleting the .gz file
			os.remove('CHECKSUMS')  # deleting the CHECKSUMS file
			gff3_to_TSSbed(newName) # this is generating a new .bed file of TSSs for the species -- saved so it can be used again
			print('\n\n{} was successfully downloaded.\n{} was successfully generated.\n\n'.format(newName, newName[:-3]+'.TSS.bed'))

##########-End-###################

#########-Start-#################
# Function: Download Fasta file from Ensembl .fa ftp site

def downloadFasta(speciesName,outputDir):
	fileList = 'wget ' + '\'' + 'ftp://ftp.ensembl.org/pub/release-98/fasta/' + speciesName + '/dna/CHECKSUMS' + '\''
	fileList_run = subprocess.run(fileList, shell = True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
	if fileList_run.returncode != 0:
		print("FAILED! Unable to locate list file. Check for existence on the Ensembl ftp site.")
		exit(2)
	else:
		with open('CHECKSUMS') as fastaList:
			for lines in fastaList:
				lines = lines.rstrip()
				if lines.endswith('dna_rm.toplevel.fa.gz'):
					temp = lines.split()
					fileGrab = temp[2]
	siteCmd = 'wget ' + '\'' + 'ftp://ftp.ensembl.org/pub/release-98/fasta/' + speciesName + '/dna/' + fileGrab + '\''
	siteCmd_run = subprocess.run(siteCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if siteCmd_run.returncode != 0:
		print("FAILED! Either your species is not available or the species did not match. Please check the Ensembl FTP site and try again")
		exit(2)
	else:
		contents=os.listdir(outputDir)
		for files in contents:
			if files.endswith('dna_rm.toplevel.fa.gz'):
				fastaGz = files
				with gzip.open(fastaGz,'rb') as f_in:
					with open(fastaGz[:-3],'wb') as f_out:
						shutil.copyfileobj(f_in, f_out)
		newName = fastaGz[:-3]
		os.remove(fastaGz)
		os.remove('CHECKSUMS')
		print('\n\n{} was successfully downloaded.'.format(newName) )

##########-End-###################

########-Start-####################
# Function: Generate user inputs for gRNA generation

def userInputs():
	print('\nWhat is your GC content lower bound (percent decimal form)?')
	gc_lower = float( input() )
	if gc_lower > 1:
		print('\nYour GC lower bound needs to be a percent decimal less than 1 and greater than or equal to 0!')
		exit(2)
	print('\nWhat is your GC content upper bound (decimal form)?')
	gc_upper = float( input() )
	if (gc_upper > 1) or (gc_upper < gc_lower):
		print('\nYour GC upper bound needs to be a decimal and greater than or equal to your lower bound!')
		exit(2)
	print('\nWhat is the minimum chromosome position for your search (integer position)?')
	minChr = int( input() )
	print('\nWhat is the maximum chromosome position for your search (integer position)?')
	maxChr = int( input() )
	if (minChr > maxChr) and (maxChr != -1):
		print('\nThe maximum chromosome position must be larger than the minimum chromosome position!')
		exit(2)
	print('\nWhat score cut-off would you like to use (range: 0-100)?')
	score = int( input() )
	if (score < 0) or (score > 100):
		print('\nThe score value must be between 0 and 100!')
		exit(2)

	inputList = [gc_lower, gc_upper, minChr, maxChr, score]
	return inputList

###########-end-####################
def main():
	parser = argparse.ArgumentParser(description="""
		The 'catPiss' program identifies a list of guide RNAs for CRISPR-based gene attenuation, promotion,\n
		or custom mode for a given list of genes and associated species. The three required inputs are a gene (or a \n
		list of genes), the species of interest (see pattern requirements below), and CRISPR mode. In addition, the \n
		user must also indicate either deNovo scoring or use of a reference scoring table. The output is a\n
		.bed file that contains a list of gRNAs for each gene and the associated primers for those gRNAs.\n 
		If you want to make deNovo gRNA (-n) from a fasta, -t and -c must be set. If you want to use your\n
		own list of gRNAs downloaded from UCSC, -f must be set to the file path.""")
	parser.add_argument('geneNames', help='The specific gene(s) of interest entered as a string or the name of the file with the list of genes.')
	parser.add_argument('speciesName', help='The name of the target species. Uses underscore for spaces (e.g., homo_sapiens)')
	parser.add_argument('targetMode' , choices = ['activation', 'interference', 'custom'], help='User determines if goal is gene activation, gene interference, or a custom mode. This argument sets the upstream and downstream boundaries relative to TSS. If "custom" is selected, then "-u" and "-d" also need to be defined')
	parser.add_argument('-o', '--output', help='Optional: supply output file directory name, otherwise write to program default', dest = 'out')
	parser.add_argument('-u', '--upstream', type=int, help='Optional: supply distance upstream of the TSS site to check for gRNAs', dest = 'upStream')
	parser.add_argument('-d', '--downstream', type=int, help='Optional: supply distance downstream of the TSS site to check for gRNAs', dest = 'downStream')
	parser.add_argument('-s', '--sortby', help="Optional: if 'de novo' is 'True' sort by is an option for sorting by score using either 'Doench2016_perc', 'Doench2016_score', 'Moreno_Matos_perc', 'Moreno_Matos_score', 'MIT_specificity'", dest = 'sortBy')
	parser.add_argument('-n', '--deNovo', help="Choice of using a de novo table from this program ('True')", dest = 'deNovoVar')
	parser.add_argument('-f', '--refTable', help="User provided reference table of scored gRNAs.", dest = 'refTableFile')
	parser.add_argument('-t', '--casType', choices = ['Cas9','Cas12A'], help="Specify type of cas: Cas9 or Cas12A. This must be entered if deNovo = 'True'", dest = 'casTypeData')
	parser.add_argument('-c', '--casFile', help="Identify file name with cas protein information detailed. This should be in the working directory. This must be entered if deNovo = 'True'", dest = 'casFileData')
	parser.add_argument('-p', '--cloning', help="Identify the cloning strategy for generation of primers. Default is pX330. Other options for cloning strategies are in 'cloningstrategy_primers.txt, which MUST be in your working directory.", dest = 'cloneDest')
	args = parser.parse_args()

	geneNames = args.geneNames
	speciesName = args.speciesName
	speciesName = speciesName.lower()
	targetModeIn = args.targetMode
	temp  = os.getcwd()
	catpath = temp + '/' + 'cat.txt' 	
	if args.deNovoVar:
		if args.casTypeData:
			casType = args.casTypeData
			casType = casType.lower()
		else:
			print("If deNovo is desired, you must define the '-t' option!")
			exit(1)
		if args.casFileData:
			casFile = args.casFileData
		else:
			print("If deNovo is desired, you must define the '-c' option!")
			exit(1)
	else: 
		if args.refTableFile:
			refTable = args.refTableFile
		else:
			print("You must provide a reference table ('-f') or choose deNovo ('-n')!")
			exit(1)
	if args.deNovoVar == False and args.refTableFile == False:
		print('You must input either a request for a deNovo analysis ("-n") or use a reference table ("-f")!')
		exit(1)
	if args.cloneDest:
		cloningStrat = args.cloneDest
	else:
		cloningStrat = 'pX330'

	currPath = os.getcwd()							# current path
	sourceDir = currPath + '/' + speciesName	# directory that contains the sources files for the species
	if args.deNovoVar:
		casFile = currPath + '/' + casFile	
		if os.path.exists(sourceDir) == 0:
			os.mkdir(sourceDir)						# if source directory for species not present, this creates it

		contents = os.listdir(sourceDir)
		TSS = False
		FASTA = False
		GRNA = False
		for files in contents:
			files = files.rstrip()
			if files.endswith('.TSS.bed'):
				TSS = True
				nameTss = files
			if files.endswith( ('.fa','.fasta') ):
				FASTA = True
				nameFasta = files
			if files.endswith('_grna_hits.bed'):
				GRNA = True
				nameGrna = files
		if (TSS == False) and (FASTA == False) and (GRNA == False):
			print("""\n\nYou are missing all of the necessary files.\nIt will take quite a bit of time to download them and\ngenerate the necessary files depending on the number of genes\n you are searching (i.e., minutes vs. days).\n\nWould you like to proceed (yes/no)?""")
			userIn1 = input()
			if userIn1 == 'yes':
				os.chdir(sourceDir)
				print('\nOK. Sit back. Relax. This is going to take a while!')
				print('\nA TSS reference file is not present. Retrieving necessary files from the Ensembl ftp download site.')
				downloadGff3(speciesName,sourceDir)
				print('\nA fasta reference file is not present. Retrieving necessary files from the Ensemble ftp download site.')
				downloadFasta(speciesName,sourceDir)
				print('\nGenerating reference guide RNA .bed file.')
				userIn = userInputs()
				denovoGuideRnaAnno(casFile, casType, nameFasta, userIn[0], userIn[1], userIn[2], userIn[3], userIn[4] )
			else:
				print('\nOk. Come back perhaps when you have more time.')
				exit(2)
		elif (TSS == True) and (FASTA == False) and (GRNA == False):
			print("""\n\nYou are missing some really large files.\nIt will take quite a bit of time to download them and\ngenerate the necessary files depending on the number of genes\n you are searching (i.e., minutes vs. days).\n\nWould you like to proceed (yes/no)?""")
			userIn1 = input()
			temp  = os.getcwd()
			catpath = temp + '/' + 'cat.txt'  
			if userIn1 == 'yes':
				os.chdir(sourceDir)
				print('\nOK. Sit back. Relax. This is going to take a while!')
				print('\nA fasta reference file is not present. Retrieving necessary files from the Ensemble ftp download site.')
				downloadFasta(speciesName,sourceDir)
				print('\nGenerating reference guide RNA .bed file.')
				userIn = userInputs()
				denovoGuideRnaAnno(casFile, casType, nameFasta,  userIn[0], userIn[1], userIn[2], userIn[3], userIn[4] )
			else:
				print('\nOk. Come back perhaps when you have more time.')
				exit(2)
		elif (TSS == True) and (FASTA == True) and (GRNA == False):
			print("""\n\nYou are missing a really big file.\nIt will take quite a bit of time to generate the necessary file depending on the number of genes\n you are searching (i.e., minutes vs. days).\n\nWould you like to proceed (yes/no)?""")
			userIn1 = input()
			if userIn1 == 'yes':
				os.chdir(sourceDir)
				print('\nOK. Sit back. Relax. This is going to take a while!')
				print('\nGenerating reference guide RNA .bed file.')
				userIn = userInputs()
				denovoGuideRnaAnno(casFile, casType, nameFasta, userIn[0], userIn[1], userIn[2], userIn[3], userIn[4] )
			else:
				print('\nOk. Come back perhaps when you have more time.')
				exit(2)
		elif (TSS == False) and (FASTA == True)  and (GRNA == False):
			print("""\n\nYou are missing a large set of the necessary files.\nIt will take quite a bit of time to download them and\ngenerate the necessary files depending on the number of genes\n you are searching (i.e., minutes vs. days).\n\nWould you like to proceed (yes/no)?""")
			userIn1 = input()
			if userIn1 == 'yes':
				os.chdir(sourceDir)
				print('\nOK. Sit back. Relax. This is going to take a while!')
				print('\nA TSS reference file is not present. Retrieving necessary files from the Ensembl ftp download site.')
				downloadGff3(speciesName,sourceDir)
				print('\nGenerating reference guide RNA .bed file.')
				userIn = userInputs()
				denovoGuideRnaAnno(casFile, casType, nameFasta, userIn[0], userIn[1], userIn[2], userIn[3], userIn[4] )
			else:
				print('Ok. Come back perhaps when you have more time.')
				exit(2)
		elif (TSS == False) and (FASTA == True) and (GRNA == True):
			print("""\n\nYou are missing a relatively small file.\nIt will take just a little bit of time (i.e., a couple of minutes) to download and\ngenerate the necessary file.\n\nWould you like to proceed (yes/no)?""")
			userIn1 = input()
			if userIn1 == 'yes':
				os.chdir(sourceDir)
				print('\nA TSS reference file is not present. Retrieving necessary files from the Ensembl ftp download site.')
				downloadGff3(speciesName,sourceDir)
			else:
				print('\nOk. Come back perhaps when you have more time.')
				exit(2)
		elif (TSS == True) and (GRNA == True):
				print('The necessary files are here.')

	if args.out:
		outputDir = currPath + '/' + outputFileDirName 
	else:
		outputDir = currPath + '/' + speciesName + '_OutputFiles'

	contents = os.listdir(sourceDir)
	for files in contents:
		if files.endswith('_grna_hits.bed'):
			nameGrna = files 

	os.chdir(currPath)
	if os.path.exists(outputDir) == 0:
		os.mkdir(outputDir)
	contents=os.listdir(sourceDir)
	for files in contents:
		if files.endswith('.TSS.bed'):
			inLocFile = sourceDir+'/' + files
			geneNamesIn = currPath + '/' + geneNames
			TSSs_of_interest_bed(inLocFile, geneNamesIn, outputDir)	# this is generating a new .bed file of TSSs specific to the desired genes
			
			overlapInput = files[:-3]+'genesofinterest.bed'
			# determining which gRNA reference table to used based on user input at command line
			if args.deNovoVar:
				gRNAs = currPath + '/' + speciesName + '/' + nameGrna
			elif args.refTableFile:
				gRNAs = refTable
			
			# initializing paramter values or using default values for 'range_overlap' function
			if targetModeIn == 'activation':
				upStrIn = 400
				downStrIn = -50
			elif targetModeIn == 'interference':
				upStrIn = 50
				downStrIn = 100
			elif targetModeIn == 'custom':
				upStrIn =  args.upStream
				downStrIn = args.downStream
			
			if args.sortBy:
				sortByIn = args.sortBy
			else:
				sortByIn = 'Doench2016_perc'
			
			if args.deNovoVar:
				deNovoIn = args.deNovoVar
			else:
				deNovoIn = False

			temp = os.getcwd()
			print('Your range for possible guides will be {} bp upstream and {} bp downstream.'.format(upStrIn, downStrIn))
			range_overlap(outputDir+'/'+overlapInput, gRNAs, outputDir+'/overLapGRnas_OUTPUT', upStrIn, downStrIn, sortByIn, deNovoIn, cloningStrat)
	 
	
	sys.exit(0)


if __name__ == '__main__':
	main()
