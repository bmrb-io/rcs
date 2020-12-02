import matplotlib.pyplot as plt
from cs_data import CS_data
from pdb_map import PDB_map
from synthesis import Synthesis
from atoms import Atom
from pdb_reader import PDB_reader
from proteins import Protein
import random
import pickle
import statistics
import time
import os

'''
cs_data = CS_data('cs_data.json')

pdb_map = PDB_map('pdb_map.json')

syn = Synthesis(cs_data.outliers_dict, pdb_map.ids_dict)
proteins_list = syn.make_proteins()

for protein in proteins_list:
    file_name = './local_proteins/' + protein.pdb_id + '_' + protein.bmrb_id + '.txt'
    print(file_name)
    picklefile = open(file_name, 'wb')
    pickle.dump(protein, picklefile)
    picklefile.close()
'''

def find_nearest_ring(atom, aromatic_rings_list):

    smallest_dist = 1e9
    smallest_ind = 0
    for i, ring in enumerate(aromatic_rings_list):
        dist = ring.find_distance(atom)
        if dist < smallest_dist:
            smallest_dist = dist
            smallest_ind = i
    if smallest_ind != 0:
        print("Damn", smallest_ind)
    return smallest_dist, smallest_ind
    

dist_list = []
shifts_list = []
angle_list = []
angle_list_down = []
angle_list_up = []

i = 0
succ = 0
fail_o = 0
fail_i = 0
npc = 0

for filename in os.listdir('./local_proteins/'):

    filename = './local_proteins/' + filename
    picklefile = open(filename, 'rb')
    protein = pickle.load(picklefile)
    protein.find_aromatic_rings()
    if len(protein.outliers_list) != 0:
        for atom in protein.outliers_list:
            if atom.position is not None and len(protein.aromatic_rings_list) > 0:
                dist, j = find_nearest_ring(atom, protein.aromatic_rings_list)
                ring = protein.aromatic_rings_list[j]
                if dist < 2e1 and len(ring.atoms_list) > 0:
                    dist_list.append(dist)
                    shifts_list.append(atom.shift_val)
                    angle = ring.find_angle(atom)
                    if angle > 90:
                        angle = 180 - angle
                    angle_list.append(angle)
                    succ += 1
                else:
                    fail_i += 1
            else:
                fail_o += 1
            if atom.position is None:
                npc += 1


print("Successfully processed:   " + str(succ))
print("Failed inner:    " + str(fail_i))
print("Failed outer:    " + str(fail_o))
print("Failed total:    " + str(fail_o + fail_i))
print("Nonetype position:    " + str(npc))

cut_off = statistics.mean(shifts_list)
for i, shift in enumerate(shifts_list):
    if shift >= cut_off:
        angle_list_up.append(angle_list[i])
    else:
        angle_list_down.append(angle_list[i]) 



fig, ax = plt.subplots()
print(len(angle_list), len(shifts_list))
#ax.scatter(angle_list, shifts_list, alpha=0.5)
ax.hist(angle_list_up, density=True, bins=20, rwidth=0.9)
ax.hist(angle_list_down, density=True, bins = 20, rwidth=0.9)
#ax.set_xlim(-1, 100)
plt.show()


