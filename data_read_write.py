import json
from noes_builder import *

def build_proteins(filename):
    proteins_dict, exceptions_k = make_proteins_from_file(filename)
    proteins_dict, exceptions_n = add_restraints(proteins_dict)
    exceptions_k.update(exceptions_n)
    exceptions_map_entries = exceptions_k
    num_entries = 0
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            num_entries += 1
    print("NUM ENTRIES:    " + str(num_entries))
    for pdb_id in exceptions_map_entries:
        for bmrb_id in exceptions_map_entries[pdb_id]:
            if pdb_id in proteins_dict:
                if bmrb_id in proteins_dict[pdb_id]:
                    del proteins_dict[pdb_id] #might be overkill to delete all pdb_id

    return proteins_dict, exceptions_map_entries

def dump_proteins(filename, proteins_dict, exceptions_map_entries):
    misaligned_list = []
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            protein = proteins_dict[pdb_id][bmrb_id]
            protein.assign_atoms_symmetrically()
            protein.prune_bad_ambiguities() #Redundant pruning
            protein.prune_missed_restraints()
            protein = protein.dump()
            proteins_dict[pdb_id][bmrb_id] = protein
    dump_dict = {
        'proteins_dict': proteins_dict,
        'exceptions_map_entries': exceptions_map_entries
    }
    with open(filename, 'w') as dumpfile:
        json.dump(dump_dict, dumpfile)

def load_proteins(filename):
    with open(filename, 'r') as dumpfile:
        dump_dict = json.load(dumpfile)
    proteins_dict = dump_dict['proteins_dict']
    exceptions_map_entries = dump_dict['exceptions_map_entries']
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            protein = Protein.load(proteins_dict[pdb_id][bmrb_id])
            proteins_dict[pdb_id][bmrb_id] = protein
    return proteins_dict, exceptions_map_entries

