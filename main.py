import os
import sys
# relative imports
from main_utils import create_dirs, check_dependencies
from assembly import get_ascensions, subset, assemble, collect_assemblies
from annotation import makeBLASTdatabase
from phylo import align, make_trees, move_trees 


# Display welcome message

message = """
mtplasmids pipeline. More info will be added 
"""
print(message)

### 0 Preparing the pipeline ###

# Function that displays missing dependencies on screen 
# with links to their webpages 
check_dependencies()

# Make working directories if they do not yet exist
# or just return them if they do
working_directory, download_directory, contigs_directory, blast_directory, annotations_directory, alignments_directory, trees_directory = create_dirs()


### 1  Downloading & de novo assembly of data ###

# Get SRA data
get_ascensions(download_directory)

# Subset downloaded read files
subset(download_directory)

# Assemble subsetted read files
assemble(download_directory)

# Collect & filter the assembly results into a single folder 
collect_assemblies(download_directory, contigs_directory)


### 2  Annotating and finding mtplasmids ###

# Make database to BLAST ORF's to
#makeBLASTdatabase(blast_directory)






### 3  Aligning & Generating trees ###

# Align fasta files obtained during annotation
#align(annotations_directory, alignments_directory)
 
# Generate maximum-likelihood trees of all generated alignments
# and move the tree to another folder to prevent cluttering
#make_trees(alignments_directory)
#move_trees(alignments_directory, trees_directory)
