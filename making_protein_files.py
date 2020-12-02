import matplotlib.pyplot as plt
from cs_data import CS_data
from pdb_map import PDB_map
from synthesis import Synthesis
from atoms import Atom
from pdb_reader import PDB_reader
from proteins import Protein
import random
import pickle
import time


#start = time.process_time()
cs_data = CS_data('cs_data.json')
#print("CSH:    " + str(time.process_time() - start))

#start = time.process_time()
pdb_map = PDB_map('pdb_map.json', cs_data.aminos_dict)
#print("MAP:    " + str(time.process_time() - start))

#start = time.process_time()
syn = Synthesis(cs_data.atoms_list, pdb_map.ids_dict)
#print("SYN:    " + str(time.process_time() - start))
start = time.process_time()
proteins_list = syn.make_proteins()

for i, protein in enumerate(proteins_list):
    file_name = './local_proteins/' + protein.pdb_id + '_' + protein.bmrb_id + '.txt'
    print(i, len(proteins_list))
    picklefile = open(file_name, 'wb')
    pickle.dump(protein, picklefile)
    picklefile.close()