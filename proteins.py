from atoms import Atom
from aromatic_rings import Aromatic_ring

class Protein:
    
    def __init__ (self, amino_index_dict, pdb_id, bmrb_id, outliers_list):
        self.pdb_id = pdb_id
        self.bmrb_id = bmrb_id
        self.amino_index_dict = amino_index_dict
        self.aromatic_rings_definition = {
            'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
        }

        self.aromatic_rings_list = []
        self.outliers_list = outliers_list

    def find_aromatic_rings(self):
        
        for atoms_list in self.amino_index_dict.values():
            amino_label = atoms_list[0].amino_label
            if amino_label in self.aromatic_rings_definition:
                ring_atoms_list = []
                atom_labels_list = self.aromatic_rings_definition[amino_label]
                for atom in atoms_list:
                    if atom.label in atom_labels_list:
                        ring_atoms_list.append(atom)
                ring = Aromatic_ring(ring_atoms_list)
                self.aromatic_rings_list.append(ring)

