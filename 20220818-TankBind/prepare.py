from Bio.PDB.PDBList import PDBList   # pip install biopython if import failure
import os
import numpy as np
import pandas as pd

pre = "./single_PDB_6hd6/"
pdir = f"{pre}/PDBs/"
pdb = '6hd6'
os.system(f"mkdir -p {pdir}")
pdbl = PDBList()
native_pdb = pdbl.retrieve_pdb_file(pdb, pdir=pdir, file_format='pdb')
