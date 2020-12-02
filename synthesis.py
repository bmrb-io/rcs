from cs_data import CS_data
from pdb_map import PDB_map
from aminos import Amino
from atoms import Atom
from pdb_reader import PDB_reader
from proteins import Protein
import random

class Synthesis:

    def __init__(self, outliers_dict, ids_dict):
        
        self.outliers_dict = outliers_dict
        self.ids_dict = ids_dict
        print(len(self.ids_dict))


    def make_proteins(self):

        proteins_list = []
        i = 0
        for bmrb_id in list(self.outliers_dict):
            print(bmrb_id)
            outliers_list = self.outliers_dict[bmrb_id]
            print(len(outliers_list))
            if bmrb_id in self.ids_dict:
                try:
                    pdbr = PDB_reader(outliers_list, bmrb_id, self.ids_dict[bmrb_id])
                    amino_index_dict = pdbr.make_amino_index_dict()
                    protein = Protein(
                        amino_index_dict, pdbr.pdb_id, pdbr.bmrb_id, outliers_list
                    )
                
                    proteins_list.append(protein)
                    i += 1
                    print(i)
                except:
                    print("Something went wrong! Oh well.")

        return proteins_list