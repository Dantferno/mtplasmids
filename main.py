import os
from assembly import get_ascensions, subset, assemble, collect_assemblies


message = """
mtplasmids pipeline. More info will be added 
"""
print(message)

# Specify the current working directories
working_directory = os.path.join('/home/jeroen/Documents/mtplasmids_pipeline')
download_directory = os.path.join(working_directory, 'sra_libraries')

# Get SRA data
#get_ascensions(download_directory)

# Subset downloaded read files
#subset(download_directory)

# Assemble subsetted read files
#assemble(download_directory)

# Collect & filter the assembly results into a single folder 
collect_assemblies(working_directory, download_directory)
