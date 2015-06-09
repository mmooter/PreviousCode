# PreviousCode

createGeneRegionsFromSim4cc.py
    - the script which creates the files for the geneRegionsFromSim4cc directory
    - can be altered for Denisovan or Neandertal by changing the file paths and names
    - output from this script:
        Files will be 4 lines.  First line will be accession number. Second line
        is the promoters in the format core|proximal|distal. Third line is exons
        separated by |.  Fourth line is introns separated by | (if there are none,
        a _ is on the fourth line).


createTranscriptsFromSim4cc.py
    - was the origial script for creating Denisovan transcripts from sim4cc output
    - can be altered for Denisovan or Neandertal by changing file paths and names

VCF2FASTA.py
    - original script for converting VCF files into fasta chromosome files
    - can be altered for Denisovan or Neandertal by changing file paths and names


sim4cc.py
    - script to run sim4cc comparison of Denisovan (or Neandertal) chromosomes
    with human transcripts from human.rna.fna
    - output will be .sim4cc.out files for each accession number with listings
    of the exons and their locations in the Denisovan (or Neandertal) chromosomes
    - can be altered for Denisovan or Neandertal by changing file paths and names

All files ending with .sh are slurm scripts for submission to the scheduling system on tusker.unl.edu
