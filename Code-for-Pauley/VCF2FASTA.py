#VCF2FASTA script v1.0
#
#Written by:
#	Benjamin Wicks
#	Jeff Kittrell
#
#Created on: 2013-03-08
#
#Modified on:
#	2013-03-12	v1.1	better documentation and attempted memory crash-recovery mechanism
#	2013-03-11	v1.0	expanded functionality and included all chromosomes
#
#This program will extract FASTA-formatted information from
#the Denisova Variant Call Format files downloaded from the
#Max Planck Institute
#(http://cdna.eva.mpg.de/denisova/VCF/hg19_1000g/)
#

import time

#constants
#TODO: Change the hard-coded filepaths
filePrefix = "AltaiNea.hg19_1000g."
fileSuffix = ".mod.vcf"
localPath = "/work/pauley/shared/projects/neanderthal/data/NeaChr/"

#fill array with all file names
r = range(1,23)
r.append('MT')
r.append('X')
r.append('Y')
#r.append('nonchrom')

#loop through all chromosome VCF files
for i in reversed(r):
	#reset counters
	numLines = 0
	numUnsures = 0
	numSkips = 0

	#if (i < 10):
	#	i = '0'+str(i)

	#open file
	fIn = open(localPath+filePrefix+str(i)+fileSuffix,"r")

	try:
		#TODO: try to open tmp file (to see where we left off)
		fIn2 = open(localPath+'FASTA/TMP/chr'+str(i)+'.tmp',"r")
		line = fIn2.readline()
		if (len(line)):
			lastChrPos = line.rstrip()
		else:
			lastChrPos = ""
		fIn2.close()
	except IOError, e:
		lastChrPos = ""

	if (len(lastChrPos) > 0):
		mode = "a"
		line1 = ""
		print("Will attempt to pick up chr "+str(i)+" after position "+lastChrPos+".")
	else:
		mode = "w"
		line1 = ">Neanderthal Chromosome "+str(i)+" built from VCF file\n"

	#ME: edited
	fOut = open(localPath+'FASTA/chr'+str(i)+'.fas',mode)
	fOut2 = open(localPath+'FASTA/TMP/chr'+str(i)+'.tmp',"w")
	fOut.write(line1)

	seqString = ""
	p1 = ""
	seqPos = 0
	startTime = time.time()

	startWriting = False

	line = fIn.readline()
	while line:
		#check for a comment (and skip over it)
		if line.startswith("#"):
			line = fIn.readline()
			continue

		#split lines at tabs
		parts = line.split("\t")
		chrPos = int(parts[1])
#		if startWriting:
#			break
			#mark down chrPos
			#fOut2 = open(localPath+'FASTA/TMP/chr'+str(i)+'.tmp',"w")
			#fOut2.write(str(chrPos))
			#fOut2.close()
#		else:
#			if (len(lastChrPos) > 0):
#				if (lastChrPos == "COMPLETED"):
#					break
#				if(int(lastChrPos) == chrPos):
#					print("Picking chr "+str(i)+" back up after position "+lastChrPos+".")
					#We are now where we left off!
#					startWriting = True
					#move along to the next line
#					line = fIn.readline()
#					continue
#				else:
					#advance and hope it gets better
#					line = fIn.readline()
#					continue
#			else:
				#We are starting from scratch on this one!
#				print("Starting chr "+str(i)+" from scratch at position "+str(chrPos)+".")
#				startWriting = True

		#read the alternate base at this position
		multiple = list(parts[3])
		alt = parts[4]

		#read the most likely base at this position
		#TODO: Use the next column to the right (ME: changed to 4 from 3)
		seqString += multiple[0]

		#put things in FASTA format
		if (len(seqString) > 80):
			p1 = seqString[:80]
			fOut.write(p1+"\n")
			p1 = ""
			seqString = seqString[80:]
		elif ((len(seqString) % 80) == 0):
			fOut.write(seqString+"\n")
			seqString = ""

		#check for gaps in the alignment (and keep a count)
		if (seqPos == 0):
			seqPos = int(parts[1])
		elif (chrPos != predicted):
			numSkips += 1

		#keep a count of how many alternate base positions there are
		if (alt != '.'):
			numUnsures += 1

		#increment the number of lines processed
		numLines += 1

		#advance line
		predicted = int(chrPos) + 1
		line = fIn.readline()

	#write remainder
	if (len(seqString)):
		fOut.write(seqString)

	#clean up
	fIn.close()	
	fOut.close()
	#ME: edited
	fOut2 = open(localPath+'FASTA/TMP/chr'+str(i)+'.tmp',"w")
	fOut2.write("COMPLETE")
	fOut2.close()

	#quick message for this file
	print "Chromosome "+str(i)+" processed in "+str(round((time.time()-startTime)/60,2))+" minutes"
	print "\tIt had "+str(numLines)+" base reads total"
	print "\tIt had "+str(numUnsures)+" base reads which had alternatives"
	print "\tIt had "+str(numSkips)+" base reads that skipped at least one genomic position\n"
