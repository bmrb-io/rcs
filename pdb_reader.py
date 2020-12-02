from atoms import Atom
import os,sys
from pdbecif.mmcif_io import CifFileReader

class PDB_reader:

    def __init__(self, outliers_list, bmrb_id, pdb_id):

        self.outliers_list = outliers_list
        self.bmrb_id = bmrb_id
        self.pdb_id = pdb_id
        print(self.pdb_id)
        self.file_path ='./cifs/{}.cif'.format(self.pdb_id)
        self.atom_site_obj = None
        self.num_models = None
        self.num_atoms = None
        if not self.check_pdb_file():
            self.get_pdb_file()
        self.make_site_obj()

    def get_pdb_file(self):
        cmd = 'wget https://files.rcsb.org/download/{}.cif -O' + self.file_path
        cmd = cmd.format(self.pdb_id)
        os.system(cmd)


    def check_pdb_file(self):
        
        file_check = os.path.isfile(self.file_path)
        
        return file_check

    def make_site_obj(self):

        cfr = CifFileReader()
        cif_obj = cfr.read(self.file_path, output='cif_wrapper')
        cif_data = cif_obj[self.pdb_id]
        self.atom_site_obj = cif_data._atom_site
        self.num_models = int(self.atom_site_obj.pdbx_PDB_model_num[-1]) 
        self.num_atoms = int(
            len(self.atom_site_obj.pdbx_PDB_model_num) / self.num_models
        )


    def make_amino_index_dict(self): #only concerned with first model for now

        atoms_list = []

        atom_label_list = self.atom_site_obj.label_atom_id[:self.num_atoms]
        pos_x_list = self.atom_site_obj.Cartn_x[:self.num_atoms]
        pos_y_list = self.atom_site_obj.Cartn_y[:self.num_atoms]
        pos_z_list = self.atom_site_obj.Cartn_z[:self.num_atoms]
        amino_label_list = self.atom_site_obj.label_comp_id[:self.num_atoms]
        amino_index_list = self.atom_site_obj.label_seq_id[:self.num_atoms]

        amino_index_dict = {}

        for i in range(self.num_atoms):
            atom = Atom()
            atom.amino_label = amino_label_list[i][:3]
            try:
                atom.amino_index = int(amino_index_list[i])
            except:
                print("Weirdness. Atom type:    " + atom_label_list[i])
                print(atom.amino_label)
            atom.label = atom_label_list[i]
            atom = self.compare_atoms(atom)
            atom.position = (
                float(pos_x_list[i]), float(pos_y_list[i]), 
                float(pos_z_list[i])
            )
            atom.index = i
            atom.is_outlier = self.compare_atoms(atom)
            atom.pdb_id = self.pdb_id
            atom.bmrb_id = self.bmrb_id

            if atom.amino_index not in amino_index_dict.keys():
                amino_index_dict[atom.amino_index] = [atom]
            else:
                amino_index_dict[atom.amino_index].append(atom)


        return amino_index_dict


    def compare_atoms(self, atom):

        for atom_outlier in self.outliers_list:
            if atom.amino_index == atom_outlier.amino_index:
                if atom.amino_label[:3] == atom_outlier.amino_label[:3]:
                    if atom.label == atom_outlier.label:
                        return atom_outlier
        
        return atom
