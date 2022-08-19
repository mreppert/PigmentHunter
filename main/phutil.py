import os,sys
import subprocess

def build_data_folders(DATADIR):
    return
    
#     if os.path.isdir(DATADIR)==0:
#         subprocess.run(['mkdir', '-p', DATADIR])

#         # Check that necessary subdirectories exist
#         subdirs = ['calc', 'exc', 'md', 'pdb', 'spec']
#         for fold in subdirs:
#             if os.path.isdir(DATADIR + "/" + fold)==0:
#                 subprocess.run(['mkdir', '-p', DATADIR + '/' + fold])

    return

def get_data_dir():
    DATADIR = os.environ['DATADIR']
    return os.path.relpath(DATADIR)

