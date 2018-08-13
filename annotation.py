import os
import subprocess
import pandas as pd
from Bio import SeqIO


# TO-DO 
# Helper function that concatenates fasta files containing DNA- & RNA-polymerase ORF's
	



# Function to generate a database of known mtplasmid DNA- & RNA-polymerase ORF's 

def makeBLASTdatabase(blast_directory):
	print('***Making ORF blast database***')
	os.chdir(blast_directory)
	
	# TO-DO
	# Incorporate helper function to concatenate DNA- & RNA-polymerase ORF sequences for the database
	
	subprocess.call('makeblastdb -in database.fsa -parse_seqids -dbtype prot -title orf_database.fsa', shell=True) 

	
# TO-DO
# Helper function that collects the right contigs and their ORF's from the original files




# Function that finds ORF's in a list of contigs and BLAST's them against a local database
# then stores the correct contigs and ORF's for further analysis

def blast_orfs(contigs_directory):
	os.chdir(contigs_directory)
	files = os.listdir(contigs_directory)
	for file in files :  
		# Store the ascension number in a variable to format output files
		ascension = 'blast{}'.format(file.strip('.')[0])	
		# Use EMBOSS's getorf to extract all ORF's from the assembled contigs
		print('*** Getting ORFS for {} ******' .format(ascension))
		subprocess.call('getorf {0} orf_{1}.fa -minsize 1000 '.format(file, ascension), shell=True) 
		# Check if the getorf output file exists to start the BLAST or if there are already results
		if os.path.exists(os.path.join(contigs_directory, 'orf_{}.fa'.format(ascension))) == True:
			# BLAST the ORFS
			print('*** BLASTING ORFS for {} ***'.format(outfile)
			blast = subprocess.check_output('blastp -db orf_database.fsa -word_size 7 -query orf_{}.fa -out blast_{} -evalue 0.1 -outfmt "10 qseqid sseqid evalue length pident" -max_target_seqs 1'.format(ascension), shell=True) 
		elif os.path.exists(os.path.join(contigs_directory, 'orf_{}.fa'.format(ascension))) == False:
			print('No ORFS were found for {}'.format(ascension))
		# TO-DO 
		# Add a check to see whether there already exists output for the contig file		



