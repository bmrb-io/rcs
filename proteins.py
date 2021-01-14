from residues import Residue
from restraints import Restraint
from atoms import Atom

class Protein:
    
    def __init__ (self, pdb_id, bmrb_id):
        self.pdb_id = pdb_id
        self.bmrb_id = bmrb_id
        self.residues_dict = {}
        self.exceptions_map_residues = {}
        self.restraints_dict = {}
        self.exceptions_map_restraints = {}

    def assign_atoms_symmetrically(self):
        exceptions_map = {}
        for restraint_id in self.restraints_dict:
            for member_id in self.restraints_dict[restraint_id]:
                restraint = self.restraints_dict[restraint_id][member_id]
                atom_amide_w_shift = self.correlate_atoms(restraint.atom_amide)
                atom_aroma_w_shift = self.correlate_atoms(restraint.atom_aroma)
                if isinstance(atom_amide_w_shift, Atom):
                    restraint.atom_amide = atom_amide_w_shift
                else:
                    exceptions_map[restraint_id] = (
                        atom_amide_w_shift
                    )
                    self.exceptions_map_restraints[restraint_id] = (
                        atom_amide_w_shift
                    )
                if isinstance(atom_aroma_w_shift, Atom):
                    restraint.atom_aroma = atom_aroma_w_shift
                else:
                    exceptions_map[restraint_id] = (
                        atom_aroma_w_shift
                    )
                    self.exceptions_map_restraints[restraint_id] = (
                        atom_aroma_w_shift
                    )

        for restraint_id in exceptions_map:
            del self.restraints_dict[restraint_id]

    def correlate_atoms(self, atom):
        if atom.res_index in self.residues_dict:
            residue = self.residues_dict[atom.res_index]
            if atom.atom_label in residue.atoms_dict:
                atom_w_shift = residue.atoms_dict[atom.atom_label]
                return atom_w_shift
            else:
                if atom.atom_label == 'H':
                    return "No such amide from k-file"
                else:
                    return atom
        else:
            return "No such residue from k-file"

    def prune_bad_ambiguities(self):
        to_prune_list = []
        for restraint_id in self.restraints_dict:
            res_index_amide_set = set()
            res_index_aroma_set = set()
            for member_id in self.restraints_dict[restraint_id]:
                restraint = self.restraints_dict[restraint_id][member_id]
                res_index_amide_set.add(restraint.atom_amide.res_index)
                res_index_aroma_set.add(restraint.atom_aroma.res_index)
            if len(res_index_amide_set) != 1 or len(res_index_aroma_set) != 1:
                to_prune_list.append(restraint_id)
                self.exceptions_map_restraints[restraint_id] = "Too many residues"
        for restraint_id in to_prune_list:
            del self.restraints_dict[restraint_id]
    
    def prune_missed_restraints(self):
        for restraint_id in self.exceptions_map_restraints:
            if restraint_id in self.restraints_dict:
                for member_id in self.restraints_dict[restraint_id]:
                    restraint = self.restraints_dict[restraint_id][member_id]
                    atom_aroma = restraint.atom_aroma
                    atom_amide = restraint.atom_amide
                del self.restraints_dict[restraint_id]
    
    def check_restraint_alignment(self):
        for restraint_id in self.restraints_dict:
            if restraint_id in self.exceptions_map_restraints:
                print(self.exceptions_map_restraints[restraint_id])
            for member_id in self.restraints_dict[restraint_id]:
                restraint = self.restraints_dict[restraint_id][member_id]
                atom_aroma = restraint.atom_aroma
                res_aroma = self.residues_dict[atom_aroma.res_index]
                if res_aroma.res_label != atom_aroma.res_label:
                    return False
                atom_amide = restraint.atom_amide
                res_amide = self.residues_dict[atom_amide.res_index]
                if res_amide.res_label != atom_amide.res_label:
                    return False
        return True
            


    def dump(self):
        dump_dict = {}
        dump_dict['pdb_id'] = self.pdb_id
        dump_dict['bmrb_id'] = self.bmrb_id
        dump_dict['residues_dict'] = {}
        for res_index in self.residues_dict:
            residue = self.residues_dict[res_index]
            dump_dict['residues_dict'][res_index] = residue.dump()
        dump_dict['exceptions_map_residues'] = self.exceptions_map_residues
        dump_dict['restraints_dict'] = {}
        for restraint_id in self.restraints_dict:
            dump_dict['restraints_dict'][restraint_id] = {}
            for member_id in self.restraints_dict[restraint_id]:
                restraint = self.restraints_dict[restraint_id][member_id]
                (
                    dump_dict['restraints_dict'][restraint_id][member_id]
                ) = restraint.dump()
        dump_dict['exceptions_map_restraints'] = self.exceptions_map_restraints
        return dump_dict

    @classmethod
    def load(cls, dump_dict):
        residues_dict = {}
        restraints_dict = {}
        pdb_id = dump_dict['pdb_id']
        bmrb_id = dump_dict['bmrb_id']
        for res_index in dump_dict['residues_dict']:
            residue = Residue.load(dump_dict['residues_dict'][res_index])
            residues_dict[res_index] = residue
        exceptions_map_residues = dump_dict['exceptions_map_residues']
        for restraint_id in dump_dict['restraints_dict']:
            restraints_dict[restraint_id] = {}
            for member_id in dump_dict['restraints_dict'][restraint_id]:
                restraint = Restraint.load(
                    dump_dict['restraints_dict'][restraint_id][member_id]
                )
                restraints_dict[restraint_id][member_id] = restraint
        exceptions_map_restraints = dump_dict['exceptions_map_restraints']
        
        
        protein = cls(pdb_id, bmrb_id)
        protein.residues_dict = residues_dict
        protein.exceptions_map_residues = exceptions_map_residues
        protein.restraints_dict = restraints_dict
        protein.exceptions_map_restraints = exceptions_map_restraints
        protein.assign_atoms_symmetrically() #this pruning is likely redundant
        protein.prune_bad_ambiguities()
        protein.prune_missed_restraints()
        return protein
