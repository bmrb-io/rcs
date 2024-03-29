from residues import Residue
from restraints import Restraint
from atoms import Atom

import traceback

class Protein:
    """
    Contain information about a protein entry in BMRB/PDB including 
    restraints and atoms. Also include methods for pruning bad restraints, etc.

    Instance variables:
    pdb_id -- the ID of the entry in PDB
    bmrb_id -- the ID of the entry in BMRB
    residues_dict -- dict containing Residue objects organized by res_index
    exceptions_map_residues -- dict of exceptions raised when creating 
        residues in k_file_reader
    restraints_dict -- dict containing Restraint objects organized by
        restraint_id and member_id
    exceptions_map_restraints -- dict of exceptions raised when creating
        restraints_dict
    pairs_dict -- dict of amide atoms and the aromatic ring protons to which
        they have a restraint

    Methods:
    __init__() -- construct an empty Protein object with only the IDs
    assign_atoms_symmetrically() -- for every restraint, assign the original 
        Atoms (with shifts)
    correlate_atoms() -- find original atoms in a Residue in residues_dict to
        replace atoms in Restraint
    prune_bad_ambiguities() -- remove from restraints_dict restraints that 
        have ambiguities to atoms from different residues
    prune_missed_restraints() -- go through exceptions_map_restraints and 
        remove any restraints still in restraints_dict
    check_restraint_alignment() -- check if the res_labels of atoms in 
        restraints_dict match those in residues_dict
    make_pairs_dict() -- make pairs_dict from restraints_dict
    prune_undefined_pairs() -- remove undefined restraints from pairs_dict
    dump() -- Create a json serializable dict containing all info in the
        Protein
    load() -- Reconstruct Protein object from a json serializable dict
    """
    
    def __init__ (self, pdb_id: str, bmrb_id: str):
        """
        Construct the empty Protein object with only the IDs.

        Keyword arguments:
        pdb_id -- the ID of the entry in PDB
        bmrb_id -- the ID of the entry in BMRB
        """

        self.pdb_id = pdb_id
        self.bmrb_id = bmrb_id
        self.residues_dict = {} # start with most of the instance variables empty dicts
        self.exceptions_map_residues = {}
        self.restraints_dict = {}
        self.exceptions_map_restraints = {}
        self.pair_geometries = {}
        self.pairs_dict = {}

    def assign_atoms_symmetrically(self):
        """
        For every restraint, assign the original Atoms (with shifts) from 
        those found in residues_dict. If not found, add an exception to 
        exceptions_map_restraints.
        """
        exceptions_map = {}
        for restraint_id in self.restraints_dict:
            for member_id in self.restraints_dict[restraint_id]:
                restraint = self.restraints_dict[restraint_id][member_id]
                # try to find the atoms from the restraint in 
                # self.residues_dict where there should be corresponding 
                # atoms _with chemical shifts_
                atom_amide_w_shift = self.correlate_atoms(restraint.atom_amide)
                atom_aroma_w_shift = self.correlate_atoms(restraint.atom_aroma)
                if isinstance(atom_amide_w_shift, Atom):
                    restraint.atom_amide = atom_amide_w_shift
                else: #probably there was no shift listed, get rid of this restraint
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

        for restraint_id in exceptions_map: #delete the restraints that didn't work
            del self.restraints_dict[restraint_id]

    def correlate_atoms(self, atom: Atom):
        """
        Find corresponding Atom object in a Residue of residues_dict. If not
        found, return an exception.

        Keyword arguments:
        atom -- the Atom object from a Restraint object to be correlated
        Returns:
        atom -- Atom object (with shift) from a Residue in residues_dict
        'No such amide from k-file' -- if the residue was found, but no atom 
            with that label
        'No such residue from k-file' -- if the residue was not found in
            residues_dict
        """
        if atom.res_index in self.residues_dict:
            residue = self.residues_dict[atom.res_index]
            if atom.atom_label in residue.atoms_dict:
                atom_w_shift = residue.atoms_dict[atom.atom_label]
                return atom_w_shift
            else:
                if atom.atom_label == 'H':
                    return "No such amide from k-file"
                else:
                    return atom #could be a pseudoatom or 5-membered TRP ring atom
        else:
            return "No such residue from k-file"

    def prune_bad_ambiguities(self):
        """
        Remove from restraints_dict restraints that have ambiguities to atoms
        from different residues.
        """
        to_prune_list = []
        for restraint_id in self.restraints_dict:
            res_index_amide_set = set()
            res_index_aroma_set = set()
            for member_id in self.restraints_dict[restraint_id]: #go through each ambiguity
                restraint = self.restraints_dict[restraint_id][member_id]
                res_index_amide_set.add(restraint.atom_amide.res_index)
                res_index_aroma_set.add(restraint.atom_aroma.res_index)
            if len(res_index_amide_set) != 1 or len(res_index_aroma_set) != 1:
                # if the restraint involves more than two res_indices, delete
                to_prune_list.append(restraint_id)
                self.exceptions_map_restraints[restraint_id] = (
                    "Too many residues"
                )
        for restraint_id in to_prune_list:
            del self.restraints_dict[restraint_id]
    
    def prune_missed_restraints(self):
        """
        Go through exceptions_map restraints. If any restraints in there are
        still in restraints_dict, remove.
        """
        for restraint_id in self.exceptions_map_restraints:
            if restraint_id in self.restraints_dict:
                for member_id in self.restraints_dict[restraint_id]:
                    restraint = self.restraints_dict[restraint_id][member_id]
                    atom_aroma = restraint.atom_aroma
                    atom_amide = restraint.atom_amide
                del self.restraints_dict[restraint_id]
    
    def check_restraint_alignment(self):
        """
        Check if the res_labels from atoms in Restraint objects of 
        restraints_dict match those in residues_dict.

        Returns:
        True -- if the labels match
        False -- if the labels do not match
        """
        not_found_count = 0
        for restraint_id in self.restraints_dict:
            for member_id in self.restraints_dict[restraint_id]:
                restraint = self.restraints_dict[restraint_id][member_id]
                atom_aroma = restraint.atom_aroma
                if atom_aroma.res_index in self.residues_dict:
                    res_aroma = self.residues_dict[atom_aroma.res_index]
                    if res_aroma.res_label != atom_aroma.res_label:
                        return False
                else:
                    self.exceptions_map_residues[atom_aroma.res_index] = (
                        "Aromatic res_index not found"
                    )
                    not_found_count += 1
                atom_amide = restraint.atom_amide
                if atom_amide.res_index in self.residues_dict:
                    res_amide = self.residues_dict[atom_amide.res_index]
                    if res_amide.res_label != atom_amide.res_label:
                        return False
                else:
                    self.exceptions_map_residues[atom_amide.res_index] = (
                        "Amide res_index not found"
                    )
                    not_found_count += 1
            if not_found_count > 10:
                return False
        return True

    def check_pair_geometries(self, cutoff=8):
        """
        Check if the mean distance (in PDB file) between any restrained 
        amide-aromatic pairs is greater than the cutoff. If not, it will be 
        disregarded in noes_builder.

        Keyword arguments:
        cutoff -- (optional, default=8) the maximum allowed distance between
            a restrained pair
        Returns:
        geom_bool -- True if all pairs are within cutoff distance, False
            otherwise
        """
        geom_bool = True
        for atom_amide in self.pairs_dict:
            geometries_amide = self.pair_geometries[atom_amide.res_index]
            for res_index_aroma in self.pairs_dict[atom_amide]:
                atoms_aroma = self.pairs_dict[atom_amide][res_index_aroma]
                labels_aroma = [atom[0].atom_label for atom in atoms_aroma]
                for atom_label in labels_aroma:
                    try:
                        if ( #otherwise, likely a pseudoatom or in 5-member TRP ring
                            atom_label in 
                            geometries_amide[res_index_aroma]
                        ):
                            pair_dist = (
                                geometries_amide[res_index_aroma][atom_label]
                            )
                            if pair_dist > cutoff:
                                geom_bool = False
                    except KeyError as err:
                        err = traceback.format_exc()
                        print(err)
                        print(labels_aroma)
                        print("HAD TO DO THIS FOR ", self.pdb_id, self.bmrb_id)
                        print(atom_amide.res_index, res_index_aroma, atom_label)
                        # there are 5 rings closer by, so this one didn't make it
                        # into an output file; if there are really just a lot of
                        # rings nearby, then it's not a problem. If there aren't, 
                        # then it will get tripped by other geometry checks.
        return geom_bool

    def make_pairs_dict(self):
        """
        Fill self.pairs_dict with amide atoms and the corresponding aromatic
        ring protons to which they have NOEs.
        """
        for restraint_id in self.restraints_dict:
            if len(self.restraints_dict[restraint_id]) == 1:
                tag = 'defi'
            elif len(self.restraints_dict[restraint_id]) > 1:
                tag = 'ambi' # so that we can ignore these later if need be
            for member_id in self.restraints_dict[restraint_id]:
                restraint = self.restraints_dict[restraint_id][member_id]
                atom_amide = restraint.atom_amide
                atom_aroma = restraint.atom_aroma
                if atom_amide not in self.pairs_dict:
                    self.pairs_dict[atom_amide] = {}
                if atom_aroma.res_index not in self.pairs_dict[atom_amide]:
                    self.pairs_dict[atom_amide][atom_aroma.res_index] = []
                self.pairs_dict[atom_amide][atom_aroma.res_index].append(
                    (atom_aroma, tag)
                )

    def prune_undefined_pairs(self):
        """
        Removes aromatic-atoms from pairs_dict that were ambiguously
        restrained to their amide atom (i.e. there were multiple members) for
        that restraint ID).
        """
        pairs_dict_new = {}
        for atom_amide in self.pairs_dict:
            pairs_dict_new[atom_amide] = {}
            for res_index_aroma in self.pairs_dict[atom_amide]:
                pairs_dict_new[atom_amide][res_index_aroma] = []
                for atom_aroma, tag in self.pairs_dict[atom_amide][res_index_aroma]:
                    if tag == 'defi':
                        pairs_dict_new[atom_amide][res_index_aroma].append(
                            (atom_aroma, tag)
                        )
        self.pairs_dict = pairs_dict_new

    def dump(self):
        """
        Create a json serializable dump_dict with all relevant info. Exclude 
        pairs_dict which can be rebuilt quickly from restraints_dict.

        Returns:
        dump_dict -- dict that can be written to json and read to reconstruct
            Protein object
        """
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
        dump_dict['pair_geometries'] = self.pair_geometries
        # just rebuild pairs_dict when in load()
        return dump_dict

    @classmethod
    def load(cls, dump_dict):
        """
        Reconstruct atom object from a json serializable dump_dict.
        
        Keyword arguments:
        dump_dict -- dict read from json and read to reconstruct Atom object
        Returns
        protein -- reconstructed Atom object
        """
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
        pair_geometries = dump_dict['pair_geometries']
        
        
        protein = cls(pdb_id, bmrb_id)
        protein.residues_dict = residues_dict
        protein.exceptions_map_residues = exceptions_map_residues
        protein.restraints_dict = restraints_dict
        protein.exceptions_map_restraints = exceptions_map_restraints
        protein.pair_geometries = pair_geometries
        protein.assign_atoms_symmetrically() #this pruning is likely redundant
        protein.prune_bad_ambiguities()
        protein.prune_missed_restraints()
        protein.make_pairs_dict() #pairs_dict was not written to dump_dict
        return protein


