class Atom:
    """
    Contain necessary information about an individual atom.

    Instance variables:
    res_index -- the index of the residue containing the atom
    res_label -- the label of the residue containing the atom
    atom_label -- the label of the atom
    cs_sigma -- the Z score of the atom

    Methods:
    __init__() -- Construct the atom and assign instance variables
    dump() -- Create a json serializable dict containing all info in the Atom
    load() -- Reconstruct Atom object from a json serializable dict
    """

    def __init__(
        self, res_index, res_label, atom_label, cs_sigma
    ):
        """
        Construct the Atom object.

        Keyword arguments:
        res_index -- the index of the residue containing the atom
        res_label -- the label of the residue containing the atom
        atom_label -- the label of the atom
        cs_sigma -- the Z score of the atom
        """
        self.res_index = res_index
        self.res_label = res_label
        self.atom_label = atom_label
        self.cs_sigma = cs_sigma

    def dump(self):
        """
        Create a json serializable dump_dict with all relevant info.

        Returns:
        dump_dict -- dict that can be written to json and read to reconstruct
            Atom object
        """
        dump_dict = {}
        dump_dict['res_index'] = self.res_index
        dump_dict['res_label'] = self.res_label
        dump_dict['atom_label'] = self.atom_label
        dump_dict['cs_sigma'] = self.cs_sigma
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
        atom_label = dump_dict['atom_label']
        cs_sigma = dump_dict['cs_sigma']
        atom = cls(res_index, res_label, atom_label, cs_sigma)
        return atom