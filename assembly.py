import os
import shutil
import subprocess
import re
# BioPython imports 
from Bio import SeqIO
from Bio.Alphabet import generic_dna


# (Optional) Function to run commands and pipe the output to the terminal 
# Only necessary if the subprocess prints to the terminal 

def run_command(command):
	p = subprocess.Popen(command, shell=True,  stdout=subprocess.PIPE)
	# Grab stdout line by line as it becomes available.  This will loop until 
	# p terminates.
	while p.poll() is None:
		l = p.stdout.readline() # This blocks until it receives a newline.
		print(l)
	# When the subprocess terminates there might be unconsumed output 
	# that still needs to be processed.
	print(p.stdout.read())


# Helper function that the deletes the original SRA database file in
# /home/user/ncbi/public/sra' after splitting to preserve drive space
  
def delete_ascension(ascension):
	ascension_folder = os.path.join('/home/{}/ncbi/public/sra'.format(getpass.getuser()))
	os.remove(os.path.join(ascension_folder,'{}.sra'.format(ascension)))
	print('.sra file of {} has been deleted!'.format(ascension))
	

# Download & split paired-end read libraries from the SRA database
# SRA Toolkit manual: https://www.ncbi.nlm.nih.gov/books/NBK242621/

def get_ascensions(download_directory):
	os.chdir(download_directory)
	# Read the ascensions from the text file
	with open('SRR_Acc_List.txt', 'r') as f:
		for line in f:
			# strip the "\n from the end of each line"
			line = line.rstrip('\n')
			# Specify the output directory
			output_path = os.path.join(download_directory, line)
			# Check if the ascension was already downloaded
			if os.path.exists(output_path) == True:
				print('Ascension {} already downloaded'.format(line))
			else:
				# Download and split the reads
				os.mkdir(output_path)
				print('Downloading & splitting: {}'.format(line))
				subprocess.call('fastq-dump -I --split-files {0} --outdir {1}'.format(line, output_path), shell=True)
				print('{} complete!'.format(line))
				# Delete the .sra file to save disk space 
				delete_ascension(line)
		print('No more ascensions to download.')

	
# Helper function to rename sequence ID's if they're improperly named during subsetting
# If the read ID's are incorrect, Abyss will throw an error and say
# that all reads are mateless

def correct_read_ids(directory, file):
	subprocess.call('sed -i "/^@{0}/ s/.2 /.1 /g" {1}'.format(directory, file), shell=True)
				
	
# Make subsets of each paired-end library of 250k reads each
# Seqtk manual: https://github.com/lh3/seqtk

def subset(download_directory):
	# Index the contents of the download directory
	paths, directories, files = next(os.walk(download_directory))
	for directory in directories:
		# Go to the sample directory and check its contents
		os.chdir(os.path.join(download_directory ,directory))
		libraries = os.listdir(os.path.join(download_directory, directory))
		# Check whether the libraries have already been subsetted
		# by checking whether the appended string "250k_" is in the filenames
		if any('250k_' in l for l in libraries):
			pass
		else: 
			# Subset the read files
			print('Subsetting {}'.format(directory))
			subprocess.call('seqtk sample -s100 {0} 250000 > 250k_{0}'.format(libraries[0]), shell=True)
			subprocess.call('seqtk sample -s100 {0} 250000 > 250k_{0}'.format(libraries[1]), shell=True)
			# Remove the old files, only the subsets are needed 
			print('Removing original files')
			os.remove(libraries[0])
			os.remove(libraries[1])
			print('250k subsets of {0} generated'.format(directory))
			libraries = os.listdir(os.path.join(download_directory, directory))
			# Correct the read id's if improperly named by seqtk
			print('Correcting read IDs')
			correct_read_ids(directory, libraries[1])
			print('Corrected reads of {}'.format(directory))
	print('All subsets have been generated')

	
# Assemble all subsetted libraries de novo using Abyss
# Abyss manual: https://github.com/bcgsc/abyss

def assemble(download_directory):
	# Index the contents of the download directory
	paths, directories, files = next(os.walk(download_directory))
	os.chdir(download_directory)
	for directory in directories:
		# Specify the path to each ascension directory & index it
		ascension_directory = os.path.join(download_directory, directory)			
		libraries = os.listdir(ascension_directory)
		# Check to see if the ascension is not yet assembled by checking
		# if only the paired-end read files are present
		if len(libraries) == 2:
			os.chdir(ascension_directory)
			subprocess.call('abyss-pe k=64 name={0} in="{1} {2}"'.format(directory, libraries[0], libraries[1]), shell=True)
			print('Ascension {} has been assembled!'.format(directory))
	print('All ascensions have been assembled!')
		

# (Optional) Helper function to extract the path of the result file containing
# the contigs from the performed assembly 	
# Only to be used if the condition in "collect_assemblies" doesn't work	
		
def get_abyss_results(directory):
	# Define a list with terms used to filter the files in the directory
	filter_list = ['scaffolds', 'contigs', 'bubbles', 'unitigs', 'indel']
	result_files = [os.path.splitext(file)[0].split('-') for file in os.listdir(directory) if file.endswith('.fa')]
	result_files[:] = [result for result in result_files if not result[1] in filter_list]
	result_files.sort()
	result_file = '{0}-{1}.fa'.format(directory, result_files[-1][1])
	return os.path.join(directory, result_file)
		

# Helper function to filter contigs in a file based on size (6000 bp < contig < 14000 bp)
# This function uses BioPython's SeqIO module to parse the file

def filter_contigs(path_to_results, filename, contigs_directory):
	# Define the path to the output file 
	path_output_file = os.path.join(contigs_directory, '{}.fa'.format(filename.split('-')[0]))
	# Check if it already exists to prevent unnecessary code execution
	if os.path.exists(path_output_file) == True:
		pass
	else:
		print('Getting results for {}'.format(filename.split('-'[0])))
		# Open the necessary files and start filtering
		with open(path_to_results, "rU") as raw_contigs, open(path_output_file, "a") as output:
			filtered_contigs = [contig for contig in raw_contigs if 6000 <= len(contig.seq) <= 14000]
			# Write to output file if contig meets criteria
			SeqIO.write(filtered_contigs, output, "fasta")					
		print('Filtered {0} contigs for ascension {1}'.format(len(filtered_contigs), filename.split('-')[0]))
	
	
# Collect & filter the assembly results in a single folder "filtered_contigs"
				
def collect_assemblies(download_directory, contigs_directory):
	# Index the contents of the download directory
	paths, directories, files = next(os.walk(download_directory))
	for directory in directories:
		# Define the necessary variables and use the "filter_contigs"
		# function to write filtered output files 
		filename = '{}-8.fa'.format(directory)
		path_to_results = os.path.join(download_directory, directory, filename)
		filter_contigs(path_to_results, filename, contigs_directory)
	print('All assembly results have been filtered!')

	
	

