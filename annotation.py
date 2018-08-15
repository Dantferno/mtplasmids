import os
import subprocess
import pandas as pd
from Bio import SeqIO


# TO-DO 
# Helper function that concatenates fasta files containing DNA- & RNA-polymerase ORF's
	
def concat_fasta(annotations_directory, blast_directory):
	dnapol = os.path.join(annotations_directory, 'dnapol.fasta')
	rnapol = os.path.join(annotations_directory, 'rnapol.fasta')
	database = os.path.join(blast_directory, 'orf_database.fasta')
	if dnapol and rnapol:
		os.chdir(annotations_directory)
		subprocess.call('cat {0} {1} > {2}'.format(dnapol, rnapol, database),shell=True)
		print('New BLAST database has been created!')
	else: 
		print('No .fasta files found! Make sure to put annotated ORFs in "annotations" as dnapol.fasta / rnapol.fasta !')

		
# Function to generate a database of known mtplasmid DNA- & RNA-polymerase ORF's 

def makeBLASTdatabase(annotations_directory, blast_directory):
	print('***Making ORF blast database***')
	concat_fasta(annotations_directory, blast_directory)
	os.chdir(blast_directory)	
	subprocess.call('makeblastdb -in orf_database.fasta -input_type fasta -dbtype prot -title orf_database.fasta', shell=True) 

	
# TO-DO
# Helper function that collects the right contigs and their ORF's from the original files




# Function that finds ORF's in a list of contigs and BLAST's them against a local database
# then stores the correct contigs and ORF's for further analysis

def blast_orfs(blast_directory, contigs_directory):
	files = os.listdir(contigs_directory)
	for file in files :  
		# Create necessary variables for the respective ascension
		ascension = file.split('.')[0]
		query = os.path.join(contigs_directory, 'orf_{}.fa'.format(ascension))
		output = os.path.join(contigs_directory, 'blast_{}.fa'.format(ascension))
		database = os.path.join(blast_directory, 'orf_database.fasta')
		# Use EMBOSS's getorf to extract all ORF's from the assembled contigs
		print('*** Getting ORFS for {} ******' .format(ascension))
		os.chdir(contigs_directory)
		subprocess.call('getorf {0} {1} -minsize 1000 -find 1 -table 4'.format(file, query), shell=True) 
		# Check if the getorf output file exists to start the BLAST or if there are already results
		if os.path.exists(os.path.join(contigs_directory, 'orf_{}.fa'.format(ascension))) == True:
			# BLAST the ORFS
			print('*** BLASTING ORFS for {} ***'.format(ascension))
			subprocess.call('blastp -db {0} -query {1} -out {2} -evalue 0.1 -outfmt "10 qseqid qlen stitle evalue length pident" -max_target_seqs 1'.format(database, query, output), shell=True) 
		elif os.path.exists(os.path.join(contigs_directory, 'orf_{}.fa'.format(ascension))) == False:
			print('No ORFS were found for {}'.format(ascension))
		# TO-DO 
		# Add a check to see whether there already exists output for the contig file		



