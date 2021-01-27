from data_read_write import *
from restraint_pair_classification import *
from noe_proportions_plotting import *
from noe_tiers import *


def make_pairs_dict(restraints_dict):
    pairs_dict = {}
    for restraint_id in restraints_dict:
        if len(restraints_dict[restraint_id]) == 1:
            tag = 'defi'
        elif len(restraints_dict[restraint_id]) > 1:
            tag = 'ambi'
        for member_id in restraints_dict[restraint_id]:
            restraint = restraints_dict[restraint_id][member_id]
            atom_amide = restraint.atom_amide
            atom_aroma = restraint.atom_aroma
            if atom_amide not in pairs_dict:
                pairs_dict[atom_amide] = {}
            if atom_aroma.res_index not in pairs_dict[atom_amide]:
                pairs_dict[atom_amide][atom_aroma.res_index] = []
            pairs_dict[atom_amide][atom_aroma.res_index].append((atom_aroma, tag))
    return pairs_dict

def prune_ambi_undef_pairs(pairs_dict):
    pairs_dict_new = {}
    for atom_amide in pairs_dict:
        pairs_dict_new[atom_amide] = {}
        for res_index in pairs_dict[atom_amide]:
            pairs_dict_new[atom_amide][res_index] = []
            for atom_aroma, tag in pairs_dict[atom_amide][res_index]:
                if  not(
                    tag == 'ambi'
                ):
                    pairs_dict_new[atom_amide][res_index].append(
                        (atom_aroma, tag)
                    )

    return pairs_dict_new

def classify_shift(cs_sigma, outlier_sigma):
    if cs_sigma <= -1 * outlier_sigma:
        return 'upfield'
    elif cs_sigma >= outlier_sigma:
        return 'downfield'
    elif abs(cs_sigma) < outlier_sigma:
        return 'non_outlier'

def results_a(dumpfile):
    """
    Gets proteins_dict and exceptions_map_entries. Prints number of amide
    protons and PDB/BMRB entries. 
    """
    proteins_dict, exceptions_map_entries = load_proteins(dumpfile)
    num_entries = 0
    num_amides = 0
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            num_entries += 1
            protein = proteins_dict[pdb_id][bmrb_id]
            for res_index in protein.residues_dict:
                residue = protein.residues_dict[res_index]
                if 'H' in residue.atoms_dict:
                    num_amides += 1
    print("A:")
    print("  NUM ENTRIES:    ", num_entries)
    print("  NUM AMIDES:    ", num_amides)
    return proteins_dict, exceptions_map_entries

def results_b1(proteins_dict, exceptions_map_entries):
    """
    Prunes empty restraints (and prints the reason). Returns the pruned proteins_dict
    """
    pruned_dict = {}
    for pdb_id in exceptions_map_entries:
        for bmrb_id in exceptions_map_entries[pdb_id]:
            exception = exceptions_map_entries[pdb_id][bmrb_id]
            if exception not in pruned_dict: ##WE Removed a del [pdb_id] here, hopefully it works
                pruned_dict[exception] = 0
            pruned_dict[exception] += 1

    pruned_dict['No amide-aromatic restraints found'] = 0
    to_del_list = []
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            protein = proteins_dict[pdb_id][bmrb_id]
            if len(protein.restraints_dict) == 0:
                pruned_dict['No amide-aromatic restraints found'] += 1
                to_del_list.append(pdb_id)
                proteins_dict[pdb_id]
    to_del_list = list(set(to_del_list))
    for pdb_id in to_del_list:
        del proteins_dict[pdb_id]
    print("B:")
    print("  I:")
    print("    " + str(pruned_dict))

    return proteins_dict

def results_b2(proteins_dict):
    """
    Takes pruned proteins dict for results_b1(). Creates pairs_dict_entries and
    returns. Counts number of entries remaining, number of amide hydrogens,
    and number of amide-aromatic pairs.
    """
    num_entries = 0
    num_pairs = 0
    pairs_dict_entries = {}
    for pdb_id in proteins_dict:
        pairs_dict_entries[pdb_id] = {}
        for bmrb_id in proteins_dict[pdb_id]:
            num_entries += 1
            protein = proteins_dict[pdb_id][bmrb_id]
            restraints_dict = protein.restraints_dict
            pairs_dict = make_pairs_dict(restraints_dict)
            for atom_amide in pairs_dict:
                for res_index in pairs_dict[atom_amide]:
                    num_pairs += 1
            pairs_dict_entries[pdb_id][bmrb_id] = pairs_dict
    
    print("  II:")
    print("    NUM ENTRIES:    ", num_entries)
    print("    NUM PAIRS:    ", num_pairs)
    return pairs_dict_entries

def results_c(proteins_dict, pairs_dict_entries):
    """
    Prunes ambiguous/undefined pairs from pairs_dict_entries and returns the
    updated dictionary. Also prints the number pruned, and the reason.
    """
    num_pairs = 0
    for pdb_id in pairs_dict_entries:
        for bmrb_id in pairs_dict_entries[pdb_id]:
            pairs_dict = pairs_dict_entries[pdb_id][bmrb_id]
            pairs_dict = prune_ambi_undef_pairs(pairs_dict)
            for atom_amide in pairs_dict:
                for res_index in pairs_dict[atom_amide]:
                    atoms_aroma = pairs_dict[atom_amide][res_index]
                    if len(atoms_aroma) != 0:
                        num_pairs += 1
            pairs_dict_entries[pdb_id][bmrb_id] = pairs_dict
    print("C:")
    print("  I:")
    print("    NUM PAIRS DEF:    ", num_pairs)
    return pairs_dict_entries

def results_d(pairs_dict_entries, outlier_sigma):

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
    for pdb_id in pairs_dict_entries:
        for bmrb_id in pairs_dict_entries[pdb_id]:
            pairs_dict = pairs_dict_entries[pdb_id][bmrb_id]
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
    return conf_tiers






    

            

#infile = 'bmrb_amide_vs_aromatic_data_his_incl.csv'
#proteins_dict, exceptions_map_entries = build_proteins(infile)
dumpfile = 'protein_dump.json' #usually protein_dump.json
#dump_proteins(dumpfile, proteins_dict, exceptions_map_entries)

outlier_sigma = 2
proteins_dict, exceptions_map_entries = results_a(dumpfile)
del_list = ['2K8V', '2KEN', '2KH9', '2KMP', '1AXH'] ###REMOVE ONE OF THE DUPLICATE PDB IDs
for pdb_id in del_list:
    bmrb_ids = list(proteins_dict[pdb_id].keys())
    bmrb_id = bmrb_ids[0]
    del proteins_dict[pdb_id][bmrb_id]
    if len(proteins_dict[pdb_id]) > 1:
        bmrb_ids = list(proteins_dict[pdb_id].keys())
        bmrb_id = bmrb_ids[0]
        del proteins_dict[pdb_id][bmrb_id]

num_duplicates = 0
for pdb_id in proteins_dict:
    if len(proteins_dict[pdb_id]) > 1:
        print("BABA BOOEY", pdb_id)
        num_duplicates += len(proteins_dict[pdb_id])
print(num_duplicates)


proteins_dict = results_b1(proteins_dict, exceptions_map_entries)
#make_proportions_plot(proteins_dict, 10, -5.5, 5.5)
pairs_dict_entries = results_b2(proteins_dict)
#make_res_prop_plot(proteins_dict, pairs_dict_entries, 10, -5.5, 5.5)
pairs_dict_entries = results_c(proteins_dict, pairs_dict_entries)

#make_res_bar_plot(pairs_dict_entries, 10, -5.5, 5.5)
conf_tiers = results_d(pairs_dict_entries, outlier_sigma)



for shift_class in conf_tiers:
    print(shift_class + ":")
    for conf in conf_tiers[shift_class]:
        print("  " + conf)
        for res_label in conf_tiers[shift_class][conf]:
            print(
                "    " + res_label + ":  " + 
                str(len(conf_tiers[shift_class][conf][res_label]))
            )

num_duplicates = 0
for pdb_id in proteins_dict:
    if len(proteins_dict[pdb_id]) > 1:
        print("BABA BOOEY", pdb_id)
        num_duplicates += len(proteins_dict[pdb_id])
print(num_duplicates)

