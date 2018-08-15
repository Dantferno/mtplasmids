import os 
import subprocess


# Function that handles all the paths to the necessary directories 

def create_dirs():
	working_directory = os.path.abspath(os.getcwd())
	
	# Specify the download directory to store SRA ascensions
	download_directory = os.path.join(working_directory, 'sra_libraries')
	if os.path.exists(download_directory) == False:
		os.mkdir(download_directory)
	
	# Specify the directory that stores the filtered contigs obtained from the assembly 
	contigs_directory = os.path.join(working_directory, 'filtered_contigs')
	if os.path.exists(contigs_directory) == False:
		os.mkdir(contigs_directory)
		
	# Specify the directory that stores blast files and databases
	blast_directory = os.path.join(working_directory, 'blast')
	if os.path.exists(blast_directory) == False:
		os.mkdir(blast_directory)
	
	# Specify the directory that stores annotations of contigs 
	annotations_directory = os.path.join(working_directory, 'annotations')
	if os.path.exists(annotations_directory) == False:
		os.mkdir(annotations_directory)
		
	# Specify the directory that stores the alignments 
	alignments_directory = os.path.join(working_directory, 'alignments')
	if os.path.exists(alignments_directory) == False:
		os.mkdir(alignments_directory)
	
	# Specify the directory that stores treefiles 
	trees_directory = os.path.join(working_directory, 'trees')
	if os.path.exists(trees_directory) == False:
		os.mkdir(trees_directory)
	
	return working_directory, download_directory, contigs_directory, blast_directory, annotations_directory, alignments_directory, trees_directory

	
# This function checks whether necessary packages have been installed
# and provides links if not 	
	
def check_dependencies():

	dependencies = os.listdir(os.path.join('/usr/bin'))
	tools = ['fastq-dump', 'seqtk', 'abyss-pe', 'muscle', 'iqtree', 'getorf']
	
	if not any(tool in tools for tool in dependencies):
		print('You are missing the following depencies. Make sure these are installed before running the pipeline!\n')
		# SRA Toolkit
		if 'fastq-dump' not in dependencies:
			subprocess.call("echo {0} '\e]8;;{1}\a{1}\e]8;;\a'".format('SRA Toolkit', 'http://ncbi.github.io/sra-tools/install_config.html'), shell=True)
		# SeqTK
		if 'seqtk' not in dependencies:
			subprocess.call("echo {0} '\e]8;;{1}\a{1}\e]8;;\a'".format('SeqTK', 'https://github.com/lh3/seqtk'), shell=True)
		# Abyss
		if 'abyss-pe' not in dependencies:
			subprocess.call("echo {0} '\e]8;;{1}\a{1}\e]8;;\a'".format('Abyss', 'https://github.com/bcgsc/abyss'), shell=True)
		# MUSCLE
		if 'muscle' not in dependencies:
			subprocess.call("echo {0} '\e]8;;{1}\a{1}\e]8;;\a'".format('MUSCLE', 'https://www.drive5.com/muscle/downloads.htm'), shell=True)
		# IQTree
		if 'iqtree' not in dependencies:
			subprocess.call("echo {0} '\e]8;;{1}\a{1}\e]8;;\a'".format('IQTree', 'http://www.iqtree.org/#download'), shell=True)
		# EMBOSS
		if 'getorf' not in dependencies:
			subprocess.call("echo {0} '\e]8;;{1}\a{1}\e]8;;\a'".format('EMBOSS', 'http://emboss.sourceforge.net/'), shell=True)
		quit()
	
	
	