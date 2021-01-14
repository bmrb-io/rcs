from atoms import Atom
from residues import Residue
from restraints import Restraint
from proteins import Protein
import os
import pynmrstar
from k_file_reader import *

def get_file(pdb_id):
    """Get restraint file from RCSB and save locally."""
    file_name = str(pdb_id).lower() + "_mr.str"
    url = "https://files.rcsb.org/download/"
    url += file_name
    cmd = 'wget ' + url + ' -O ./noe_rcsb_strs/' + file_name
    os.system(cmd)

def get_star_restraints(pdb_id):
    """
    Download restraint file if not available locally, and check for issues with
    the file.

    Keyword arguments:
    pdb_id -- string of the PDB ID
    Returns:
    'No restraint file' -- if no file can be found from RCSB
    'Bad restraint file' -- if file downloaded cannot be parsed by pynmrstar
    'No restraints in file' -- if there are no Gen_dist_constrain loops in file
    'Unacceptable restraint loop type' -- if an unexpected distance constraint
        subtype is found
    restraint_loops_list -- list of restraint loops from restraint file
    """
    filename = f"{pdb_id.lower()}_mr.str"
    filepath = os.path.join('noe_rcsb_strs', filename)
    if not os.path.isfile(filepath):
        get_file(pdb_id)
    try:
        entry = pynmrstar.Entry.from_file(filepath)
    except AttributeError:
        return "No restraint file"
    except ValueError:
        return "Bad restraint file"
    restraint_loops_list = entry.get_loops_by_category("Gen_dist_constraint")
    if len(restraint_loops_list) == 0:
        return "No restraints in file"
    if check_noe_loops(entry):
        return restraint_loops_list
    else:
        return 'Unexpected restraint loop subtype'

def check_noe_loops(entry):
    """
    Loop through Constraint_file loop in restraint file and check if only
    expected restraint loop subtypes with type distance are in restraint file.

    Keyword arguments:
    entry -- a pynmrstar entry
    Returns:
    True -- if no unexpected subtypes found
    False -- if any unexpected subtypes found
    """
    info_loop = entry.get_loops_by_category("Constraint_file")[0]
    info_list = info_loop.get_tag(
        ["Constraint_type", "Constraint_subtype"]
    )
    num_noe_loops = 0
    good_subtypes = [
        'NOE', 'general distance', 'hydrogen bond', 'disulfide bond', 'PRE'
    ]
    for info in info_list:
        if info[0] == 'distance':
            subtype = info[1]
            if subtype not in good_subtypes:
                return False
    return True

def check_amide(atom_1, atom_2):
    """
    Check if either of atom_1 or atom_2 is an amide hydrogen; if so, return
    that atom.

    Keyword arguments:
    atom_1 -- the first Atom object built from the restraint in file
    atom_2 -- the second Atom object built from the restraint in file
    Returns:
    bool_amide -- True if one of the atoms is an amide, False otherwise
    atom_amide -- None if neither of the atoms is an amide, otherwise it is the
       atom that was found to be an amide
    """
    atom_amide = None
    bool_amide = False
    atoms_list = [atom_1, atom_2]
    for atom in atoms_list:
        if atom.atom_label == 'H':
            bool_amide = True
            atom_amide = atom
    return bool_amide, atom_amide

def check_aromatic(atom):
    """
    Check if an atom is an aromatic ring proton.

    Keyword arguments:
    atom -- the Atom object to be checked
    Returns:
    bool_aroma -- True if the atom is an aromatic ring proton, False otherwise
    atom_aroma -- The Atom object
    """
    bool_aroma = False
    atom_aroma = None
    aromatics_dict = {
        "PHE": ["HD1", "HD2", "QD", "HE1", "HE2", "QE", "HZ"],
        "TYR": ["HD1", "HD2", "QD", "HE1", "HE2", "QE", "HH"],
        "TRP": ["HD1", "HE1", "HE3", "QE", "HZ2", "HZ3", "QZ", "HH2"],
        "HIS": ["HD1", "HD2", "QD", "HE1", "HE2", "QE"]
    }
    if atom.res_label in aromatics_dict:
        if atom.atom_label in aromatics_dict[atom.res_label]:
            bool_aroma = True
            atom_aroma = atom
    return bool_aroma, atom_aroma

def make_restraint(restraint_entry):
    """
    Read line from restraint file and build restraint if it is amide-aromatic.

    Keyword arguments:
    restraint_entry -- str of a line of restraint data in the restraint file
    Returns:
    'No amide atom or atoms in same residue' -- self-explanatory exception
    'No aromatic ring proton' -- if one atom is amide but other is not aromatic
        ring proton
    restraint -- Restraint object with an amide atom and an aromatic ring 
        proton
    restraint_id -- restraint ID as found in restraint file
    member_id -- member ID as found in restraint file
    """

    restraint_id = restraint_entry[0]
    member_id = restraint_entry[1]
    logic_code = restraint_entry[2]

    res_index_1 = restraint_entry[3]
    res_label_1 = restraint_entry[4]
    atom_label_1 = restraint_entry[5]
    atom_1 = Atom(res_index_1, res_label_1, atom_label_1, None)  

    res_index_2 = restraint_entry[6]
    res_label_2 = restraint_entry[7]
    atom_label_2 = restraint_entry[8]
    atom_2 = Atom(res_index_2, res_label_2, atom_label_2, None)

    bool_amide, atom_amide = check_amide(atom_1, atom_2)
    if bool_amide and atom_1.res_index != atom_2.res_index:
        if atom_1 == atom_amide:
            bool_aroma, atom_aroma = check_aromatic(atom_2)
        elif atom_2 == atom_amide:
            bool_aroma, atom_aroma = check_aromatic(atom_1)
        else:
            print(bool_amide, atom_amide)
        if bool_aroma:
            restraint = Restraint(
                atom_amide, atom_aroma
            )
        else:
            return "No aromatic ring proton", restraint_id, member_id
    else:
        return "No amide atom or atoms in same residue", restraint_id, member_id

    return restraint, restraint_id, member_id

def make_restraints_dict(pdb_id):
    """
    Create dictionary of all restraints for a PDB entry.

    Keyword arguments:
    pdb_id -- string of PDB ID for entry
    Returns
    restraints_dict -- dict of all restraints for entry by restraint_id and 
        member_id
    exceptions_map -- dict (by restraint_id) of exceptions raised by
        make_restraint()
    """
    restraints_dict = {}
    exceptions_map = {}
    restraint_loops_list = get_star_restraints(pdb_id)
    if not isinstance(restraint_loops_list, list):
        return restraint_loops_list, None
    for i, restraint_loop in enumerate(restraint_loops_list):
        restraints_list = restraint_loop.get_tag(
            [
                "ID", "Member_ID", "Member_logic_code",
                "Comp_index_ID_1", "Comp_ID_1", "Atom_ID_1",
                "Comp_index_ID_2", "Comp_ID_2", "Atom_ID_2",
                "Distance_val", "Distance_lower_bound_val", 
                "Distance_upper_bound_val"
            ]
        )
        for restraint_entry in restraints_list:
            restraint, restraint_id, member_id = make_restraint(restraint_entry)
            restraint_id = str(i) + ',' + restraint_id #to tell which restraint loop it is from!
            if isinstance(restraint, Restraint):
                if restraint_id not in restraints_dict:
                    restraints_dict[restraint_id] = {}
                restraints_dict[restraint_id][member_id] = restraint
            else:
                exceptions_map[restraint_id] = restraint
            
    
    return restraints_dict, exceptions_map


def add_restraints(proteins_dict): 
    """
    Add restraints_dict to all possible proteins in proteins_dict; call Protein
    methods to correlate atoms in restraints with atoms from BMRB; prune 
    various bad restraints.

    Keyword arguments:
    proteins_dict -- dict of all proteins by pdb_id and bmrb_id
    Reutrns
    proteins_dict -- now with restraints added to those possible
    exceptions_map_entry -- dict (by pdb_id and bmrb_id) of exceptions raised
        when building proteins and their restraints_dicts
    """
    exceptions_map_entries = {}
    for pdb_id in proteins_dict:
        print(pdb_id)
        for bmrb_id in proteins_dict[pdb_id]:
            protein = proteins_dict[pdb_id][bmrb_id]
            restraints_dict, exceptions_map_restraints = make_restraints_dict(
                pdb_id
            )
            if not isinstance(restraints_dict, dict):
                if pdb_id not in exceptions_map_entries:
                    exceptions_map_entries[pdb_id] = {}
                exceptions_map_entries[pdb_id][bmrb_id] = restraints_dict
            else:
                protein.restraints_dict = restraints_dict
                protein.assign_atoms_symmetrically()
                protein.prune_bad_ambiguities()
                protein.prune_missed_restraints()
                if protein.check_restraint_alignment():
                    protein.exceptions_map_restraints = exceptions_map_restraints
                    proteins_dict[pdb_id][bmrb_id] = protein
                else:
                    if pdb_id not in exceptions_map_entries:
                        exceptions_map_entries[pdb_id] = {}
                    exceptions_map_entries[pdb_id][bmrb_id] = (
                        "Misaligned restraint indices"
                    )

    return proteins_dict, exceptions_map_entries

