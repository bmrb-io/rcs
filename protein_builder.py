from k_file_maker import *
from k_file_reader import *
from noes_builder import *
from ids_map import *

def build_protein(pdb_id, bmrb_id):
    """
    Create a Protein object with amides, aromatic protons, chemical shifts,
    and NOE restraints if possible. Otherwise return an appropriate exception.

    Keyword arguments:
    pdb_id -- str of the PDB ID
    bmrb_id -- str of the BMRB ID
    Returns
    protein -- either a Protein object or an exception raised by 
        make_protein_from_file() or add_restraints()
    """
    ring_current_object = RingCurrentEffect(pdb_id, bmrb_id)
    k_file_path = ring_current_object.calculate_ring_current_effects(
        pdb_id, bmrb_id
    )
    protein = make_protein_from_file(k_file_path)
    if isinstance(protein, Protein):
        protein = add_restraints(protein)
    return protein


def build_all_proteins(map_filename):
    """
    Create a Protein object for all viable entries, and return an exception 
    otherwise. Add successful proteins to proteins_dict and exceptions to
    exceptions_map.

    Keyword arguments:
    map_filename -- path to file containing map of bmrb_ids to pdb_ids
    Returns:
    proteins_dict -- dict of all successfully created proteins
    exceptions_map -- dict of reasons for failure to create proteins

    """
    id_map = IDs_map(map_filename)
    proteins_dict = {}
    exceptions_map = {}
    for bmrb_id, pdb_id in id_map.ids_list[:50]:
        protein = build_protein(pdb_id, bmrb_id)
        if isinstance(protein, Protein):
            print("Success!", pdb_id, bmrb_id)
            if pdb_id not in proteins_dict:
                proteins_dict[pdb_id] = {}
            proteins_dict[pdb_id][bmrb_id] = protein
        else:
            print("Failure!", protein, pdb_id, bmrb_id)
            if pdb_id not in exceptions_map:
                exceptions_map[pdb_id] = {}
            exceptions_map[pdb_id][bmrb_id] = protein 