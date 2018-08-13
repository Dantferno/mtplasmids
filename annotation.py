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
	# Incorporate helper function to concatenate ORF sequences for the database
	
	subprocess.call('makeblastdb -in database.fsa -parse_seqids -dbtype prot -title orf_database.fsa', shell=True) 


	
# Function that collects all found ORF's in a file for easy lookup 

def collect_orfs(file):
	# Create the dataframe 
	df = pd.Dataframe()
	# Append all ORF's to the dataframe 
	print('*** Collecting ORFS for {} ***'.format(file))
	with open('ORFtmp.fa', 'rU') as f:
		for record in SeqIO.parse(f, "fasta"): 
			df = df.append({'Description': record.description, 'Sequence': record.seq}, ignore_index=True) #Append sequence and description of ORF to DataFrame
		print('Number of ORFs found : {}' .format(df.shape[0]))
		# Store the dataframe in a .csv file
		df.to_csv('{}_orfs.csv'.format(file.strip('.')[0]))

	
# TO-DO
# Helper function that collects the right contigs and their ORF's from the original files




# Function that finds ORF's in a list of contigs and BLAST's them against a local database
# then stores the correct contigs and ORF's for further analysis

def blast_orfs(contigs_directory):
	os.chdir(contigs_directory)
	files = os.listdir(contigs_directory)
	for file in files :                 
		# Use EMBOSS's getorf to extract all ORF's from the assembled contigs
		print('***Getting ORFS for {}******' .format(file))
		subprocess.call('getorf {} ORFtmp.fa -minsize 1000 '.format(file), shell=True) 
		# Check if the getorf outout file exists to start the BLAST 
		if os.path.exists(os.path.join(contigs_directory, 'ORFtmp.fa')) == True:
			# Collect all the ORF's in .csv file 
			collect_orfs(file)
			# BLAST the ORFS
			print('*** BLASTING ORFS for {} ***'.format(file.strip('.')[0]))
			outfile = 'blast{}'.format(file.strip('.')[0])
			blast = subprocess.check_output('blastp -db orf_database.fsa -word_size 7 -query ORFtmp.fa -out {} -evalue 0.1 -outfmt "10 qseqid sseqid evalue length pident" -max_target_seqs 1'.format(outfile), shell=True) 
			os.remove('ORFtmp.fa')



