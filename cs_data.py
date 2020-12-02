import requests
import json
from aminos import Amino
from atoms import Atom
import time

class CS_data:
    
    def __init__(self, file_name, make_file = False):

        self.file_name = file_name
        self.aminos_id_list = [
            'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS',
            'LEU', 'MET', 'ASN', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP',
            'TYR'
        ] #PRO excluded
        self.aminos_dict = {}
        self.outliers_dict = {}
        if make_file:
            self.make_cs_data_json()
        start = time.process_time()
        cs_data_json = self.get_cs_data_json()
        print("JSON:    " + str(time.process_time() - start))

        start = time.process_time()
        self.make_aminos_dict(cs_data_json)
        self.make_outliers_dict()
        print("DICT:    " + str(time.process_time() - start))
    
    def make_aminos_dict(self, cs_data_json):

        self.init_aminos_dict()
        for atom_entry in cs_data_json['data']:
            #print(atom_entry)
            atom = Atom()
            atom.bmrb_id = atom_entry[0]
            atom.amino_label = atom_entry[3][:3]
            atom.amino_index = atom_entry[2]
            atom.shift_val = float(atom_entry[6]) 
            atom.label = 'H'
            if atom.amino_label in self.aminos_id_list:
                amino = self.aminos_dict[atom.amino_label]
                amino.atoms_list.append(atom)

        for amino in self.aminos_dict.values():
            amino.make_gaussian_fit()
            amino.find_outliers(2, 4)


    def init_aminos_dict(self):

        for amino_label in self.aminos_id_list:
            amino = Amino(amino_label)
            self.aminos_dict[amino_label] = amino
        
    def get_cs_data_json(self):
        
        with open(self.file_name) as f:
            cs_data_json = json.load(f)

        return cs_data_json

    def make_cs_data_json(self):
    
        cs_params_dict = {}
        cs_params_dict['atom_id'] = 'H'
        cs_url = "r."
        r = requests.get(cs_url, cs_params_dict).json()
        with open(self.file_name, 'w') as f:
            json.dump(r, f)

    def make_outliers_dict(self):

        for amino in self.aminos_dict.values():
            for atom in amino.atoms_list:
                if atom.bmrb_id not in self.outliers_dict:
                    self.outliers_dict[atom.bmrb_id] = []
                if atom.is_outlier:
                    self.outliers_dict[atom.bmrb_id].append(atom)

