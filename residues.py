from atoms import Atom
class Residue:
    """
    Contain necessary information about a residue, most importantly atoms_dict.

    Instance variables:
    res_index -- the index of the residue within the entry
    res_label -- the label of the residue within the entry
    atoms_dict -- dict of the atoms within the residue, organized by atom label

    Methods:
    __init__() -- Construct the Residue and assign instance variables
    dump() -- Create a json serializable dict containing all info in the Residue
    load() -- Reconstruct Residue object from a json serializable dict
    """
    def __init__(self, res_index, res_label, atoms_dict):
        """
        Construct the Residue object.

        Keyword arguments:
        res_index -- the index of the residue
        res_label -- the label of the residue
        atoms_dict -- dict containing Atom objects organized by atom label
        """
        self.res_index = res_index
        self.res_label = res_label
        self.atoms_dict = atoms_dict
    
    def dump(self):
        """
        Create a json serializable dump_dict with all relevant info.

        Returns:
        dump_dict -- dict that can be written to json and read to reconstruct
            Residue object
        """
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
        """
        Reconstruct atom object from a json serializable dump_dict.
        
        Keyword arguments:
        dump_dict -- dict read from json and read to reconstruct Atom object
        Returns
        atom -- reconstructed Atom object
        """
        res_index = dump_dict['res_index']
        res_label = dump_dict['res_label']
        atoms_dict = {}
        for atom_label in dump_dict['atoms_dict']:
            atom = Atom.load(dump_dict['atoms_dict'][atom_label])
            atoms_dict[atom_label] = atom
        residue = cls(res_index, res_label, atoms_dict)
        return residue
