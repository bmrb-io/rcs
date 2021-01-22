from ids_map import *
from protein_builder import *
from proteins import *
import requests

def get_proteins_dict(entries_dict, build_anyway=False):
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
    proteins_dict = {}
    exceptions_map = {}
    for pdb_id in entries_dict:
        for bmrb_id in entries_dict[pdb_id]:
            print(pdb_id)
            protein = get_protein(pdb_id, bmrb_id, build_anyway)
            if isinstance(protein, Protein):
                if pdb_id not in proteins_dict:
                    proteins_dict[pdb_id] = {}
                proteins_dict[pdb_id][bmrb_id] = protein
            else:
                if pdb_id not in exceptions_map:
                    exceptions_map[pdb_id] = {}
                exceptions_map[pdb_id][bmrb_id] = protein
    return proteins_dict, exceptions_map

def get_all_entries():
    """Return a dict of all corresponding PDB and BMRB IDs."""
    url = "http://api.bmrb.io/v2/mappings/bmrb/pdb"
    r = requests.get(url).json()
    entries_dict = {}
    i = 0
    for ids_dict in r:
        i += 1
        if i <= 10:
            bmrb_id = ids_dict['bmrb_id']
            pdb_ids = ids_dict['pdb_ids']
            for pdb_id in pdb_ids:
                if pdb_id not in entries_dict:
                    entries_dict[pdb_id] = []
                entries_dict[pdb_id].append(bmrb_id)
        else:
            break
    return entries_dict