import os
from assembly import get_ascensions, subset, assemble, collect_assemblies
from phylo import align, make_trees, move_trees 

message = """
mtplasmids pipeline. More info will be added 
"""
print(message)

# Specify the current working directories
working_directory = os.path.join('/home/jeroen/Documents/mtplasmids_pipeline')
download_directory = os.path.join(working_directory, 'sra_libraries')
annotations_directory = os.path.join(working_directory, 'annotations')
alignments_directory = os.path.join(working_directory, 'alignments')
trees_directory = os.path.join(working_directory, 'trees')

### 1  Downloading & de novo assembly of data ###

# Get SRA data
#get_ascensions(download_directory)

# Subset downloaded read files
#subset(download_directory)

# Assemble subsetted read files
#assemble(download_directory)

# Collect & filter the assembly results into a single folder 
#collect_assemblies(working_directory, download_directory)


### 2  Annotating the downloaded data ###


### 3  Aligning & Generating trees ###

# Align fasta files obtained during annotation
#align(annotations_directory, alignments_directory)
 
# Generate maximum-likelihood trees of all generated alignments
# and move the tree to another folder to prevent cluttering
#make_trees(alignments_directory)
move_trees(alignments_directory, trees_directory)
