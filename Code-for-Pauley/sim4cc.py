#!/usr/bin/python
import sys
sys.path.append('/work/pauley/shared/programs/tools/biopython-1.56/')   
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import os
import subprocess

fIn = open('/work/pauley/shared/projects/Neandertal/human.rna.fna','r')
fIn2Prefix = '/work/pauley/shared/projects/denisovan_lab/data/humanGenomeFASTAChr/'
for record in SeqIO.parse(fIn, 'fasta'):
        parts = record.id.split('|')
        theID = parts[1]
        accNum = parts[3]
	#Check if it has already been sim4cc'd
        fExisting = '/work/pauley/shared/projects/Neandertal/rawSim4cc/'+str(accNum)+'.sim4cc.out'
        if os.path.isfile(fExisting):
                statinfo = os.stat(fExisting)
                if (statinfo.st_size > 0):
                        print 'Already sim4cc\'d '+accNum+' - skipping.'
                        continue
        fIn2Name = fIn2Prefix+accNum+'.chr.txt'
	fComparePath = '/work/pauley/shared/projects/Neandertal/AltaiNeaGenome.fna'
        if os.path.isfile(fIn2Name):
                statinfo = os.stat(fIn2Name)
                if (statinfo.st_size > 0):
        		print "Trying "+accNum+"..."
			fIn2 = open(fIn2Name,'r')
			chrNum = fIn2.readline().strip()
			fIn2.close()
			if (chrNum == 'Unknown'):
				#chrNum = 'MT'
				print accNum+' on unknown chromosome - skipping.'
				#Remove temporary file before exiting iteration
				continue
			fComparePath = '/work/pauley/shared/projects/Neandertal/chromosomesFromVCF/chr'+chrNum+'.fas'
			cmd = ['/work/pauley/shared/programs/src/sim4cc.2010-11-22/sim4cc', '/work/pauley/shared/projects/denisovan_lab/data/humanGenomeFASTA/'+accNum+'.fasta', '/work/pauley/shared/projects/Neandertal/chromosomesFromVCF/chr'+chrNum+'.fas', 'A=0', 'R=2']
			return_code = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = return_code.communicate()
			fOut = open('/work/pauley/shared/projects/Neandertal/rawSim4cc/'+accNum+'.sim4cc.out','w')
			fOut.write(out)
			fOut.close()
			if (len(err)):
				fOut2 = open('/work/pauley/shared/projects/Neandertal/rawSim4cc/'+accNum+'.sim4cc.err','w')
				fOut2.write(err)
				fOut2.close()
			print 'sim4cc complete for '+accNum+'! Output saved'
		else:
			print 'No chromosome information for '+accNum+' - skipping.'

	print 'Compare Path for '+accNum+': '+fComparePath
