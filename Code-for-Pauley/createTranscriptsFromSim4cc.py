import sys
sys.path.append('/work/pauley/shared/programs/tools/biopython-1.56/')
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import os

#Get input files and set paths to files
inPath = "/work/pauley/shared/projects/Denisova/rawSim4cc/"
inputFile = open("/work/pauley/shared/projects/denisovan_lab/data/human.rna.fna", "r")
outPath = "/work/pauley/shared/projects/Denisova/transcriptsFromSim4cc/"

#Iterate through all records in human.rna.fna, using accession numbers to make exon files
for record in SeqIO.parse(inputFile, 'fasta'):
	i = 0
	exon = []
	complement = 0
	#Retrieve accession number
        parts = record.id.split('|')
        theID = parts[1]
        accNum = parts[3]
	inExisting = '/work/pauley/shared/projects/Denisova/rawSim4cc/'+str(accNum)+'.sim4cc.out'
	fExisting = '/work/pauley/shared/projects/Denisova/transcriptsFromSim4cc/'+str(accNum)+'.sim4cc.exon'
	if os.path.isfile(fExisting):
		statinfo = os.stat(fExisting)
		if (statinfo.st_size > 0):
			print "Skipping " +str(accNum) + ". . ."
			continue
	if not(os.path.isfile(inExisting)):
		continue
	#Use accession to open input and output file
	inFile2 = open((str(inPath)+accNum+".sim4cc.out"), "r")
	outFile = open(str(outPath)+str(accNum)+".sim4cc.exon", "w")
	print "Beginning parse of " + str(accNum)
	#Read and save lines of input file
	blank = inFile2.readline()
	seq1 = inFile2.readline()
	seq2 = inFile2.readline()
	line = inFile2.readline()
	while (line):
		if (line == "\n"):
			line = inFile2.readline()
                        continue
		if (line.startswith("(complement)")):
			line = inFile2.readline()
			complement = 1
			continue
		else:
			line2 = line.split("  ")
			#if ('==' in str(line2[2])):
			#	line = inFile2.readline()
			#	continue
			exon.append(line2[1])
			line = inFile2.readline()
			i = i + 1
	if (len(exon) == 0):
		print "No exon matches found. Skipping this record."
		continue
	print "Input file read, analyzing now"
	#Get file name for Denisovan chromosome and open
	chrNum = seq2.split(" ")
	chrFile = open(str(chrNum[2]),"r")
        fastaHeader = chrFile.readline()
        line = chrFile.readline()
        chrLine = ""
        while (line):
                line1 = line.rstrip()
                chrLine = chrLine + line1
                line = chrFile.readline()
	#Iterate through exons and write to file
	j = 0
	seq = ""
	while (j < len(exon)):
		#Format exon ends and save
		chrOneLine = chrLine
		range = exon[j].replace( "(", "")
		range1 = range.replace( ")", "")
		ends = range1.split("-")
		startPos = ends[0]
		counter = int(startPos) - 1
		endPos = ends[1]
		#Open file
		#chrFile = open(str(chrNum[2]),"r")
		#fastaHeader = chrFile.readline()
		#line = chrFile.readline()
		#chrOneLine = ""
		#while (line):
		#        line1 = line.rstrip()
        	#	chrOneLine = chrOneLine + line1
        	#	line = chrFile.readline()
		while (counter < int(endPos)):
			seq = seq + str(chrOneLine[counter])
			counter = counter + 1
		j = j + 1
	print "Exons read and added to sequence.  Writing to file now"
	#Get reverse complement if complement=1
	if (complement == 1):
		revComp = Seq(str(seq))
		seq = revComp.reverse_complement()
	#Write results to file and close files
	outFile.write(">Denisovan transcript " + str(accNum) +  ", from " + str(chrNum[2]) + "\n")
	numLines = (len(seq)/80) + 1
        x = 0
        while (x < numLines):
        	if (len(seq)>=80):
                	p1 = seq[:80]
                	outFile.write(str(p1) + "\n")
                        p1 = ""
                        seq = seq[80:]
                else:
                	outFile.write(str(seq))
                x = x + 1
	outFile.close()
	chrFile.close()
	inFile2.close()
	print "Parse of " + str(accNum) + "COMPLETE"
inputFile.close()
