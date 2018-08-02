import os
import shutil


# Download & split the paired-end read libraries from the SRA database
# SRA Toolkit manual: https://www.ncbi.nlm.nih.gov/books/NBK242621/

def get_ascensions(download_directory):
	os.chdir(download_directory)
	with open('sra_list.txt', 'r') as f:
		for line in f:
			output_path = os.path.join(download_directory, line)
			print('Downloading & splitting: {}'.format(line))
			os.system('fastq-dump -I --split-files {0} --outdir {1}'.format(line, output_path))
			print('{} complete!'.format(line))
