from atoms import Atom
from typing import Dict

class Restraint:
    
    def __init__(self, atom_amide: Atom, atom_aroma: Atom):
        self.atom_amide = atom_amide
        self.atom_aroma = atom_aroma

    def dump(self) -> Dict:

        dump_dict = {}
        dump_dict['atom_amide'] = self.atom_amide.dump()
        dump_dict['atom_aroma'] = self.atom_aroma.dump()
        return dump_dict

    @classmethod
    def load(cls, dump_dict):

        atom_amide = Atom.load(dump_dict['atom_amide'])
        atom_aroma = Atom.load(dump_dict['atom_aroma'])
        restraint = cls(atom_amide, atom_aroma)

        return restraint