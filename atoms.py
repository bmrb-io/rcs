class Atom:
    """
    Contain necessary information about an individual atom.

    Instance variables:
    
    """

    def __init__(
        self, res_index, res_label, atom_label, cs_sigma
    ):
        self.res_index = res_index
        self.res_label = res_label
        self.atom_label = atom_label
        self.cs_sigma = cs_sigma

    def dump(self):

        dump_dict = {}
        dump_dict['res_index'] = self.res_index
        dump_dict['res_label'] = self.res_label
        dump_dict['atom_label'] = self.atom_label
        dump_dict['cs_sigma'] = self.cs_sigma
        return dump_dict
    
    @classmethod
    def load(cls, dump_dict):

        res_index = dump_dict['res_index']
        res_label = dump_dict['res_label']
        atom_label = dump_dict['atom_label']
        cs_sigma = dump_dict['cs_sigma']
        atom = cls(res_index, res_label, atom_label, cs_sigma)
        return atom