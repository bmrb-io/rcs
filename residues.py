from atoms import Atom
class Residue:

    def __init__(self, res_index, res_label, atoms_dict):
        self.res_index = res_index
        self.res_label = res_label
        self.atoms_dict = atoms_dict
    
    def dump(self):
        dump_dict = {}
        dump_dict['res_index'] = self.res_index
        dump_dict['res_label'] = self.res_label
        dump_dict['atoms_dict'] = {}
        for atom_label in self.atoms_dict:
            atom = self.atoms_dict[atom_label]
            dump_dict['atoms_dict'][atom_label] = atom.dump()
        return dump_dict

    @classmethod
    def load(cls, dump_dict):
        res_index = dump_dict['res_index']
        res_label = dump_dict['res_label']
        atoms_dict = {}
        for atom_label in dump_dict['atoms_dict']:
            atom = Atom.load(dump_dict['atoms_dict'][atom_label])
            atoms_dict[atom_label] = atom
        residue = cls(res_index, res_label, atoms_dict)
        return residue
