import os
import subprocess


# Function to make alignements out of fasta files 
# obtained during annotation of the assembly 

def align(annotations_directory, alignments_directory):
	# Index the contents of the annotations directory
	paths, directories, files = next(os.walk(annotations_directory))
	os.chdir(annotations_directory)
	# Align each file using MUSCLE 
	for file in files:
		print('Aligning {}'.format(file))
		output_path = os.path.join(alignments_directory, file.split('.')[0])
		subprocess.call('muscle -in {0} -out {1}.fasta'.format(file, output_path), shell=True)
	print('All {} files have been aligned!'.format(len(files)))
		
	
# Function to generate maximum-likelihood trees for 
# the previously made alignments

def make_trees(alignments_directory):
	# Index the contents of the annotations directory
	paths, directories, files = next(os.walk(alignments_directory))
	os.chdir(alignments_directory)
	# Generate a maximum-likelihood tree for each alignment
	for file in files:
		print('generating tree for {}'.format(file))
		if '_aa' in file: 
			subprocess.call('iqtree -s {0} -pre tree_{1} -m VT+F+I+G4 -alrt 1000 -bb 1000 -nt AUTO'.format(file, file.split('.')[0]), shell=True)
		elif '_dna' in file:
			subprocess.call('iqtree -s {0} -pre tree_{1} -alrt 1000 -bb 1000 -nt AUTO'.format(file, file.split('.')[0]), shell=True)
	print('All trees have been generated!')
	
	
# Function to relocate the tree files 
# to prevent cluttering of the alignments folder	
	
def move_trees(alignments_directory, trees_directory):
	# Index the contents of the annotations directory
	paths, directories, files = next(os.walk(alignments_directory))
	os.chdir(alignments_directory)
	# Find all alignment files in the directory
	alignments = [file for file in files if '.fasta' in file]
	for alignment in alignments:
		# Find all the treefiles belonging to each alignment
		treefiles = [file for file in files if 'tree_{}'.format(alignment.split('.')[0]) in file.split('.')[0]]
		# If tere are treefiles for the alignment, creata directory and move them
		if len(treefiles) > 0:
			output_path = os.path.join(trees_directory, '{}'.format(alignment.split('.')[0]))
			if not output_path:
				os.mkdir(output_path)
			for file in treefiles:
				os.rename(file, os.path.join(output_path, file))
			print('Treefiles for {} have been moved'.format(alignment.split('.')[0]))
		else:
			print('No treefiles found for {}'.format(alignment.split('.')[0]))
	print('All treefiles have been moved')