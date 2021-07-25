import json

from en_masse import *
from noe_proportions_plotting import *
from typing import Union, Dict

def classify_shift(cs_sigma: float, outlier_sigma: Union[int, float]) -> str:
    """
    Classify the chemical shift depending on the number of standard deviations
    that define an outlier.

    Keyword arguments:
    cs_sigma -- the Z-score of an atom
    outlier_sigma -- the number of standard deviations from mean that make an
        atom's chemical shift an outlier
    Returns:
    'upfield' -- if the shift is upfield by more than outlier_sigma standard 
        deviations
    'downfield' -- if the shift is downfield by more than outlier_sigma 
        standard deviations
    'non_outlier' -- if the shift is within outlier_sigma standard deviations
        of the mean
    """
    if cs_sigma <= -1 * outlier_sigma:
        return 'upfield'
    elif cs_sigma >= outlier_sigma:
        return 'downfield'
    elif abs(cs_sigma) < outlier_sigma:
        return 'non_outlier'

def results_a(
    proteins_dict: Dict[str, Dict[str, Protein]], 
    exceptions_map_entries: Dict[str, Dict[str, str]]
):
    """
    The first stage of printing out the results. Prints
    exceptions_map_entries in more readable format; prints the number of
    entries that had restraints successfully added; prints the number of
    amide-aromatic pairs with at least one NOE between them.

    Keyword arguments:
    proteins_dict -- dict organized by pdb_id and bmrb_id of all succesfully
        created Protein instances
    exceptions_map_entries -- dict organized by pdb_id and bmrb_id of all 
        exceptions raised from failures to add restraints
    """
    print("A:")
    expected_exceptions = [
        'No matching atoms found', 'Too many entities/assemblies', 
        'No distance restraints in file', 'Unreleased structure', 'No pairs found',
        'No restraint file', 'Empty restraint file', 'Too many restraints',
        'Misaligned restraint indices', 'BMRB entry only exists in NMR-STAR 2.0',
        'No aromatic residues found', 'Misformatted restraint file',
        'Unacceptable distances between restrained pairs',
        'Restraint file not in reboxitory',
        'mmCIF file not in reboxitory',
        'STR file not in reboxitory',
        'Permission denied error (this should no longer occur)'
    ]
    exceptions_by_reason = {}
    for pdb_id in exceptions_map_entries:
        for bmrb_id in exceptions_map_entries[pdb_id]:
            reason = exceptions_map_entries[pdb_id][bmrb_id]
            if reason not in exceptions_by_reason: 
                #print(f"{reason}: {bmrb_id}, {pdb_id}")
                exceptions_by_reason[reason] = 0
            exceptions_by_reason[reason] += 1 # Add to counter of exceptions with this reason
            if reason not in expected_exceptions: #if there was an unanticipated exception
                print(
                    f"UNEXPECTED EXCEPTION IN {pdb_id}, {bmrb_id}: {reason}"
                )
                pass
    for reason in exceptions_by_reason: 
        num = exceptions_by_reason[reason]
        print("  ", reason, ":", num)
    num_entries = 0
    num_pairs = 0
    num_amides = 0
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            num_entries += 1
            protein = proteins_dict[pdb_id][bmrb_id]
            pairs_dict = protein.pairs_dict
            for atom_amide in pairs_dict:
                for res_index_aroma in pairs_dict[atom_amide]:
                    num_pairs += 1


    print("  ", "NUM ENTRIES WITH USABLE RESTRAINTS:    ", num_entries)
    print("  ", "NUM AMIDE-AROMATIC PAIRS WITH A RESTRAINT:    ", num_pairs)



def results_b(proteins_dict: Dict[str, Dict[str, Protein]]):
    """
    The second stage of printing out the results. Prunes undefined pairs from
    each protein.pairs_dict. Prints the number of amide-aromatic pairs with at
    least one defined restraint between them.

    Keyword arguments:
    proteins_dict -- dict organized by pdb_id and bmrb_id of all successfully
        created Protein instances
    """
    print("B:")
    num_pairs = 0
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            protein = proteins_dict[pdb_id][bmrb_id]
            protein.prune_undefined_pairs() # We don't want ambiguous pairs here. 
            pairs_dict = protein.pairs_dict
            for atom_amide in pairs_dict:
                for res_index_aroma in pairs_dict[atom_amide]:
                    atoms_aroma = pairs_dict[atom_amide][res_index_aroma]
                    if len(atoms_aroma) != 0: # Some will be empty because they only had ambiguous restraints
                        num_pairs += 1
    print(
        "  ", "NUM AMIDE_AROMATIC PAIRS WITH A DEFINED RESTRAINT:    ", num_pairs
    )

    print_restraint_exceptions(proteins_dict)

def basic_write(proteins_dict, pairs_to_write, filename):
    ids = []
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            ids.append([pdb_id, bmrb_id])
    
    with open(f'{filename}_ids.json', 'w') as df:
        json.dump(ids, df)
    
    with open(f'{filename}_pairs.json', 'w') as df:
        json.dump(pairs_to_write, df)

def print_restraint_exceptions(proteins_dict):
    exc_by_reason = {}

    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            protein = proteins_dict[pdb_id][bmrb_id]
            exceptions_map = protein.exceptions_map_restraints
            for restraint_id in exceptions_map:
                reason = exceptions_map[restraint_id]
                if reason not in exc_by_reason:
                    exc_by_reason[reason] = 0
                exc_by_reason[reason] += 1
    
    print('RESTRAINT EXCEPTIONS:')
    for reason in exc_by_reason:
        print(f"  {reason}:    {exc_by_reason[reason]}")

def print_result_stages(
    outlier_sigma: Union[int, float], build_anyway: bool = False, 
    make_plots: bool = False, only_recoord: bool = False,
):
    """
    Generates dict of all proteins with BMRB and PDB entries, analyzes 
    restraint data, and prints out major results.

    Keyword arguments:
    outlier_sigma -- the number of standard deviations from mean that make an
        atom's chemical shift an outlier
    build_anyway -- default False; if True, proteins will be built even if a 
        json build file is present in ./proteins
    """

    proteins_dict, exceptions_map_entries = get_proteins_dict_multi(
        build_anyway
    )
    
    if make_plots:
        make_all_plots(proteins_dict, 11, -5.5, 5.5)
        make_num_restraints_plot(proteins_dict, outlier_sigma)
    results_a(proteins_dict, exceptions_map_entries) 

    results_b(proteins_dict)


print_result_stages(2, make_plots=False)
