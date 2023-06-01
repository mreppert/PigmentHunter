import os

# Build data folder structure if necessary. 
def build_data_folders(DATADIR):
    
    # Check if DATADIR exists
    if os.path.isdir(DATADIR)==0:
        # If not, create it. 
        subprocess.run(['mkdir', '-p', DATADIR])

        # Now create any necessary subdirectories
        subdirs = ['calc', 'exc', 'md', 'pdb', 'spec']
        for fold in subdirs:
            if os.path.isdir(DATADIR + "/" + fold)==0:
                subprocess.run(['mkdir', '-p', DATADIR + '/' + fold])
    return