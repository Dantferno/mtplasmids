import os 


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
	
	return working_directory, download_directory, contigs_directory, annotations_directory, alignments_directory, trees_directory
	

