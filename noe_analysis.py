from en_mass import *
from restraint_pair_classification import *
from noe_proportions_plotting import *
from noe_tiers import *

def classify_shift(cs_sigma, outlier_sigma):
    if cs_sigma <= -1 * outlier_sigma:
        return 'upfield'
    elif cs_sigma >= outlier_sigma:
        return 'downfield'
    elif abs(cs_sigma) < outlier_sigma:
        return 'non_outlier'

def results_a(proteins_dict, exceptions_map_entries):
    print("A:")
    exceptions_by_reason = {}
    for pdb_id in exceptions_map_entries:
        for bmrb_id in exceptions_map_entries[pdb_id]:
            reason = exceptions_map_entries[pdb_id][bmrb_id]
            if reason not in exceptions_by_reason:
                exceptions_by_reason[reason] = 0
            exceptions_by_reason[reason] += 1
    for reason in exceptions_by_reason:
        num = exceptions_by_reason[reason]
        print("  ", reason, ":", num)
    num_entries = 0
    num_pairs = 0
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

def results_b(proteins_dict):
    """
    """
    print("B:")
    num_pairs = 0
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            protein = proteins_dict[pdb_id][bmrb_id]
            protein.prune_ambi_undef_pairs()
            pairs_dict = protein.pairs_dict
            for atom_amide in pairs_dict:
                for res_index_aroma in pairs_dict[atom_amide]:
                    atoms_aroma = pairs_dict[atom_amide][res_index_aroma]
                    if len(atoms_aroma) != 0:
                        num_pairs += 1
    print(
        "  ", "NUM AMIDE_AROMATIC PAIRS WITH A DEFINED RESTRAINT:    ", num_pairs
    )

def results_c(proteins_dict, outlier_sigma):

    conf_tiers = {
        'upfield': {
            'high': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []},
            'intermediate': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []},
            'low': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []}, 
            'none': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []}
        },
        'downfield': {
            'high': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []},
            'intermediate': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []},
            'low': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []}, 
            'none': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []}
        },
        'non_outlier': {
            'high': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []},
            'intermediate': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []},
            'low': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []}, 
            'none': {'HIS': [], 'TRP': [], 'PHE': [], 'TYR': []}
        },
    }
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            protein = proteins_dict[pdb_id][bmrb_id]
            protein.prune_ambi_undef_pairs()
            pairs_dict = protein.pairs_dict
            for atom_amide in pairs_dict:
                for res_index in pairs_dict[atom_amide]:
                    atoms_aroma = pairs_dict[atom_amide][res_index]
                    atoms_aroma = [atom[0] for atom in atoms_aroma] #removing tags
                    if len(atoms_aroma) > 0:
                        conf_tier, label_aroma = get_confidence_tier(
                            atoms_aroma
                        )
                        shift_class = classify_shift(
                            atom_amide.cs_sigma, outlier_sigma
                        )
                        conf_tiers[shift_class][conf_tier][label_aroma].append(
                            [pdb_id, bmrb_id, atom_amide, atoms_aroma]
                        )


    print("D:")
    for shift_class in conf_tiers:
        print(f"  {shift_class}:")
        for conf_tier in conf_tiers[shift_class]:
            print(f"    {conf_tier}:")
            for res_label in conf_tiers[shift_class][conf_tier]:
                num = conf_tiers[shift_class][conf_tier][res_label]
                print(f"      {res_label}:  {num}")


def print_result_stages(
    outlier_sigma, build_anyway=False
    ):
    entries_dict = get_all_entries()
    proteins_dict, exceptions_map_entries = get_proteins_dict(
        entries_dict, build_anyway
    )
    results_a(proteins_dict, exceptions_map_entries)
    results_b(proteins_dict)
    results_c(proteins_dict, outlier_sigma)
