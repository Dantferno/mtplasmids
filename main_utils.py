import os
import shutil
import subprocess


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


# Download & split paired-end read libraries from the SRA database
# SRA Toolkit manual: https://www.ncbi.nlm.nih.gov/books/NBK242621/

def get_ascensions(download_directory):
	os.chdir(download_directory)
	# Read the ascensions from the text file
	with open('SRR_Acc_List.txt', 'r') as f:
		for line in f:
			# Specify the output directory
			output_path = os.path.join(download_directory, line)
			# Check if the ascension was already downloaded
			if os.path.exists(output_path) == True:
				print('Ascension {} already downloaded'.format(line))
			else:
				# Download and split the reads
				os.mkdir(output_path)
				print('Downloading & splitting: {}'.format(line))
				subprocess.call('fastq-dump -I --split-files {0} --outdir {1}'.format(line, output_path))
				print('{} complete!'.format(line))
		print('No more ascensions to download.')

			
# Make subsets of each paired-end library of 250k reads each
# Seqtk manual: https://github.com/lh3/seqtk

def subset(download_directory):
	# Index the contents of the working directory
	paths, directories, files = next(os.walk(download_directory))
	for directory in directories:
		# Go to the sample directory and check its contents
		os.chdir(os.path.join(download_directory ,directory))
		libraries = os.listdir(os.path.join(download_directory, directory))
		# Check whether the libraries have already been subsetted
		# by checking whether the appended string "250k_" is in the filenames
		if any('250k_' in l for l in libraries):
			print('Ascension {} has already been subsetted'.format(directory))
		else: 
			# Subset the read files
			subprocess.call('seqtk sample -s100 {0} 250000 > 250k_{0}'.format(libraries[0]))
			subprocess.call('seqtk sample -s100 {0} 250000 > 250k_{0}'.format(libraries[1]))
			# Remove the old files, only the subsets are needed 
			os.remove(libraries[0])
			os.remove(libraries[1])
			print('250k subsets of sample {0} generated'.format(directory))
	print('All subsets have been generated')

	
# Assemble all subsetted libraries de novo using Abyss
# Abyss manual: https://github.com/bcgsc/abyss

def assemble(download_directory):
	# Index the contents of the working directory
	paths, directories, files = next(os.walk(download_directory))
	os.chdir(download_directory)
	for directory in directories:
		# Specify the path to each ascension directory 
		ascension_directory = os.path.join(download_directory, directory)
		# Check whether the ascension has already been assembled
		if os.path.exists(os.path.join(ascension_directory, '{}-stats.csv'.format(directory) == True:
			print('Ascension {} has already been assembled!'.format(directory))
		else:
			libraries = os.listdir(ascension_directory)
			# Second check to see if assembly is necessary 
			# and run, if so
			if libraries == 2:
				os.chdir(ascension_directory)
				subprocess.call('abyss-pe k=64 name={0} in="{1} {2}"'.format(directory, libraries[0], libraries[1]))
				print('Ascension {} has been assembled!'.format(directory))
	print('All ascensions have been assembled!')
		
	

