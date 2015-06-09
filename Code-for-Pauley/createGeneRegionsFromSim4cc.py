import sys
sys.path.append('/work/pauley/shared/programs/tools/biopython-1.56/')
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import os

#Get input files and set paths to files
inPath = "/work/pauley/shared/projects/Denisova/rawSim4cc/"
inputFile = open("/work/pauley/shared/projects/denisovan_lab/data/human.rna.fna", "r")
#outPath = "/work/pauley/shared/projects/Denisova/geneRegionsFromSim4cc/"
outPath = "/work/pauley/shared/projects/Scripts/tempGeneRegions/"

#Iterate through all records in human.rna.fna, using accession numbers to make output files
for record in SeqIO.parse(inputFile, 'fasta'):
        i = 0
        exon = []
	intron = []
        complement = 0
        #Retrieve accession number
        parts = record.id.split('|')
        accNum = parts[3]
        inExisting = '/work/pauley/shared/projects/Denisova/rawSim4cc/'+str(accNum)+'.sim4cc.out'
        fExisting = '/work/pauley/shared/projects/Denisova/geneRegionsFromSim4cc/'+str(accNum)+'.sim4cc.final'
	    if os.path.isfile(fExisting):
                statinfo = os.stat(fExisting)
                if (statinfo.st_size > 0):
                        print "Skipping " +str(accNum) + ". . ."
                        continue
        if not(os.path.isfile(inExisting)):
                continue
        #Use accession to open input and output file
        inFile2 = open((str(inPath)+accNum+".sim4cc.out"), "r")
        outFile = open(str(outPath)+str(accNum)+".sim4cc.final", "w")
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
                        exon.append(line2[1])
                        line = inFile2.readline()
                        i = i + 1
        if (len(exon) == 0):
                print "No exon matches found. Skipping this record."
                continue
        print "Input file read, analyzing now"
        #Get file name for Denisovan chromosome and open
        chrNum = seq2.split(" ")
	chrNum = chrNum[2]
	chrom = ""
	if (chrNum == "/work/pauley/shared/databases/VCF/T_hg19_1000g/FASTA/chrX.fas" or chrNum == "/work/pauley/shared/databases/VCF/T_hg19_1000g/FASTA/chrY.fas"):
		chrom = chrNum[-8:]
	else:
		chrom = chrNum[-9:]
	outFile.write(">Denisovan " + str(accNum) +  ", from " + str(chrom) + "\n")
        chrFile = open("/work/pauley/shared/projects/Denisova/chromosomesFromVCF/"+str(chrom),"r")
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
        nonCompStart = 0
        CompStart = 0
        while (j < len(exon)):
                #Format exon ends and save
                range = exon[j].replace( "(", "")
                range1 = range.replace( ")", "")
                ends = range1.split("-")
                startPos = ends[0]
                endPos = ends[1]
		if (j == 0):
                        intron.append(endPos)
                if (j != 0):
                        intron.append(startPos)
                        intron.append(endPos)
                if (complement == 0 and j == 0):
                        nonCompStart = startPos
                        break
                if (complement == 1 and j == (len(exon)-1)):
                        CompStart = endPos
                j = j + 1
        dist = 0
	distal = ""
	proximal = ""
	core = ""
        #Retrieve promoters when not complement
        if (complement == 0):
                dist = int(nonCompStart) - 3000
                #Get promoters and write to files
                counter = dist - 1
                while (counter < int(nonCompStart) - 1):
                    chrOneLine = chrLine
                    seq = seq + str(chrOneLine[counter])
                    counter = counter + 1
		distal = seq
		proximal = seq[-1000:]
		core = seq[-300:]

        #Get reverse complement if complement=1
        if (complement == 1):
                dist = int(CompStart) + 3000
                counter = int(CompStart)
                while (counter < dist):
                    chrOneLine = chrLine
                    seq = seq + str(chrOneLine[counter])
                    counter = counter + 1
		revComp = Seq(str(seq))
		seq = revComp.reverse_complement()
		distal = seq
		proximal = seq[-1000:]
		core = seq[-300:]

	#Next print promoters to file.  Then gather exon sequences and print them to file.  Following that, gather intron sequences and print them to file
	#Format as follows: >accnum\ncorepromoter|proximalpromoter|distalpromoter\nexon|exon|....\nintron|intron|...
	#If no introns (single exon), then place _ on line
	outFile.write(str(core)+"|"+str(proximal)+"|"+str(distal)+"\n")
	print "PROMOTERS WRITTEN"
	j = 0
        exonSeqs = ""
	seq = ""
        while (j < len(exon)):
                #Format exon ends and save
                chrOneLine = chrLine
		seq = ""
                range = exon[j].replace( "(", "")
                range1 = range.replace( ")", "")
                ends = range1.split("-")
                startPos = ends[0]
                counter = int(startPos) - 1
                endPos = ends[1]
                while (counter < int(endPos)):
                        seq = seq + str(chrOneLine[counter])
                        counter = counter + 1
		#Get reverse complement if complement=1
		if (complement == 1):
			revComp = Seq(str(seq))
			seq = revComp.reverse_complement()
		if (j != (len(exon)- 1)):
			exonSeqs = exonSeqs + str(seq) + "|"
		if (j == (len(exon)-1)):
			exonSeqs = exonSeqs + str(seq) + "\n"
                j = j + 1
	outFile.write(exonSeqs)
	print "EXONS WRITTEN"
	intronSeqs = ""
	seq = ""
	if (len(exon) == 1):
		intronSeqs = "_"
	else:
		h = 0
		intron.sort()
		while (h < (len(intron) - 1)):
			counter = int(intron[h])
			seq = ""
                	chrOneLine = chrLine
                	while (counter < ((int(intron[h+1]))-1)):
                	    seq = seq + str(chrOneLine[counter])
                	    counter = counter + 1
			#Get reverse complement if complement=1
                	if (complement == 1):
                        	revComp = Seq(str(seq))
                        	seq = revComp.reverse_complement()
			if (h != (len(intron)-2)):
				intronSeqs = intronSeqs + str(seq) + "|"
			if (h == (len(intron)-2)):
				intronSeqs = intronSeqs + str(seq) + "\n"
			h = h + 2
	outFile.write(intronSeqs)
	print "INTRONS WRITTEN"
	outFile.close()
        chrFile.close()
        inFile2.close()
        print "Parse of " + str(accNum) + " COMPLETE"
inputFile.close()
