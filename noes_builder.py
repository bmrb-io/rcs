from atoms import Atom
from residues import Residue
from restraints import Restraint
from proteins import Protein
from k_file_reader import *
import os
import pynmrstar
import requests
from typing import Union, List, Tuple, Dict


def get_star_restraints(pdb_id: str) -> Union[List, str]: #find out type in List
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

    filepath = os.path.join(
        '/reboxitory', '2021', '06', 'PDB', 'data', 'structures', 'all', 
        'nmr_restraints_v2', f"{pdb_id.lower()}_mr.str.gz"
    )
    try:
        entry = pynmrstar.Entry.from_file(filepath)
    except AttributeError as err:
        return "Empty restraint file"
    except pynmrstar.exceptions.ParsingError:
        return "No restraint file"
    except ValueError:
        return "Misformatted restraint file"
    restraint_loops_list = entry.get_loops_by_category("Gen_dist_constraint")
    if len(restraint_loops_list) == 0:
        return "No distance restraints in file" 
    loops_check = check_noe_loops(entry)
    if loops_check == 'All clear':
        return restraint_loops_list
    else:
        return loops_check

def check_noe_loops(entry) -> str: #find out type of entry
    """
    Loop through Constraint_file loop in restraint file and check if: only
    expected restraint loop subtypes with type distance are in restraint file; 
    no loops have more than 3500 restraints, which may indicate issues with 
    the file.

    Keyword arguments:
    entry -- a pynmrstar entry
    Returns:
    'All clear' -- if no unexpected subtypes found and loops are not too long
    'Too many restraints' -- if any distance restraint loops are too long
    'Undexpected restraint_loop_subtype' -- if any unexpected subtypes found
    """
    info_loop = entry.get_loops_by_category("Constraint_file")[0]
    info_list = info_loop.get_tag(
        ["Constraint_type", "Constraint_subtype", "Constraint_number"]
    ) 
    good_subtypes = [
        'NOE', 'general distance', 'hydrogen bond', 'disulfide bond', 'PRE'
    ] # the expected subtypes
    for info in info_list:
        if info[0] == 'distance':
            subtype = info[1]
            if int(info[2]) > 3500:
                return "Too many restraints"
            if subtype not in good_subtypes:
                return "Unexpected restraint_loop_subtype"
    return "All clear"

def check_amide(atom_1: Atom, atom_2: Atom) -> Tuple[bool, Union[Atom, None]]:
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

def check_aromatic(atom: Atom) -> Tuple[bool, Union[Atom, None]]:
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

def check_dist(dist_val: str, dist_lower: str, dist_upper: str) -> str:
    """
    Check that distance values and bounds of restraint are acceptable.

    Keyword arguments:
    dist_val -- the distance reported for the restraint
    dist_lower -- the lower bound of the distance
    dist_upper -- the upper bound of the distance
    Returns:
    'Only lower bound reported' -- if dist_val and dist_upper are not included 
        for the restraint
    'No upper, dist_val too high' -- if no dist_upper, and dist_val is too high
    'Upper too high' -- if dist_upper is too high, regardless of others
    'All clear' -- if dist_upper is not too high or (if dist_upper is not 
        reported) if dist_val is not too high
    """
    if dist_upper == '.':
        if dist_val == '.':
            return "Only lower bound reported" # This is indicative that something is wrong
        else:
            if float(dist_val) <= 5:
                return "All clear"
            else:
                return "No upper, dist_val too high" # I have not seen this triggered
    else:
        if float(dist_upper) <= 6:
            return "All clear"
        else:
            return "Upper too high" # Again, indicative that something is wrong, maybe not an NOE

def make_restraint(restraint_entry) -> Tuple[Union[Restraint, str], str, str]: #find out type of restraint_entry
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
    # info for cataloguing restraints
    restraint_id = restraint_entry[0]
    member_id = restraint_entry[1]
    logic_code = restraint_entry[2]
    # info for the first atom
    res_index_1 = restraint_entry[3]
    res_label_1 = restraint_entry[4]
    atom_label_1 = restraint_entry[5]
    atom_1 = Atom(res_index_1, res_label_1, atom_label_1, None, None)  
    # info for the second atom
    res_index_2 = restraint_entry[6]
    res_label_2 = restraint_entry[7]
    atom_label_2 = restraint_entry[8]
    atom_2 = Atom(res_index_2, res_label_2, atom_label_2, None, None)
    #info for check_dist
    dist_val = restraint_entry[9]
    dist_lower = restraint_entry[10]
    dist_upper = restraint_entry[11]

    bool_amide, atom_amide = check_amide(atom_1, atom_2)
    if bool_amide and atom_1.res_index != atom_2.res_index:
        if atom_1 == atom_amide:
            bool_aroma, atom_aroma = check_aromatic(atom_2)
        elif atom_2 == atom_amide:
            bool_aroma, atom_aroma = check_aromatic(atom_1)
        if bool_aroma:
            dist_check = check_dist(dist_val, dist_lower, dist_upper)
            if dist_check == "All clear":
                restraint = Restraint(
                    atom_amide, atom_aroma
                )
            else:
                return dist_check, restraint_id, member_id
        else:
            return "No aromatic ring proton", restraint_id, member_id
    else:
        return "No amide atom or atoms in same residue", restraint_id, member_id

    return restraint, restraint_id, member_id

def make_restraints_dict(
    pdb_id: str
) -> Tuple[Dict[str, Dict[str, Restraint]], Dict[str, Dict[str, str]]]:
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
    if not isinstance(restraint_loops_list, list): # exception triggered
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
            else: # exception triggered
                exceptions_map[restraint_id] = restraint
    return restraints_dict, exceptions_map

def add_restraints(protein: Protein) -> Union[Protein, str]:
    """
    Add restraints dict to protein.

    Keyword arguments:
    protein -- Protein object with amide and aromatic ring protons included.
    Returns:
    restraints_dict -- dictionary of Restraint objects organized by restraint
        ID and member ID
    'No restraint file' -- if no file can be found from RCSB
    'Bad restraint file' -- if file downloaded cannot be parsed by pynmrstar
    'No restraints in file' -- if there are no Gen_dist_constrain loops in file
    'Unacceptable restraint loop type' -- if an unexpected distance constraint
        subtype is found
    'No pairs found' -- if no amide-aromatic restraints were found
    'Misaligned restraint indices' -- if the residues have different indices in
        PDB file and restraint file
    
    """
    pdb_id = protein.pdb_id
    restraints_dict, exceptions_map_restraints = make_restraints_dict(pdb_id)
    if not isinstance(restraints_dict, dict):
        #returns exception thrown by make_restraints_dict()
        return restraints_dict 
    else:
        protein.restraints_dict = restraints_dict
        protein.exceptions_map_restraints = exceptions_map_restraints
        if protein.check_restraint_alignment():
            protein.assign_atoms_symmetrically()
            protein.prune_bad_ambiguities()
            protein.prune_missed_restraints()
            if len(restraints_dict) == 0:
                # successfully built restraints_dict, but no acceptable amide-aromatic restraints
                return "No pairs found"
            protein.make_pairs_dict()
            if protein.check_pair_geometries():
                return protein
            else:
                return "Unacceptable distances between restrained pairs"
        else:
            return "Misaligned restraint indices"

