import os
from main_utils import get_ascensions, subset, assemble


message = """
mtplasmids pipeline. More info will be added 
"""
print(message)

# Specify the current working directory
working_directory = os.path.join('/home/jeroen/Documents/mtplasmids_pipeline/sra_libraries')

# Get SRA data
get_ascensions(working_directory)

# Subset downloaded read files
subset(working_directory)

# Assemble subsetted read files
assemble(working_directory)
