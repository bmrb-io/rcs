import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd

def make_binning_info(num_bins, cs_min, cs_max):

    bin_edges = list(np.linspace(cs_min, cs_max, num_bins))
    bin_midpoints = [
        (bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(bin_edges) - 1)
    ]
    
    return bin_edges, bin_midpoints

def bin_atoms(atoms_list, bin_edges, bin_midpoints):

    cs_list = []
    for atom in atoms_list:
        cs_sigma = float(atom.cs_sigma)
        if abs(cs_sigma) <= bin_edges[-1]:
            cs_list.append(cs_sigma)
    inds = np.digitize(cs_list, bin_edges)
    binned_list = [0 for i in bin_midpoints]
    for bin_ind in inds:
        bin_ind = bin_ind - 1
        binned_list[bin_ind] += 1
    
    return binned_list

def sort_atoms_by_restraint(protein):
    restraintless_set = set()
    restraintful_set = set()
    for restraint_id in protein.restraints_dict:
        member_id0 = list(protein.restraints_dict[restraint_id].keys())[0]
        restraint = protein.restraints_dict[restraint_id][member_id0]
        amide_atom = restraint.atom_amide
        restraintful_set.add(amide_atom)
    for res_index in protein.residues_dict:
        residue = protein.residues_dict[res_index]
        if 'H' in residue.atoms_dict:
            atom_amide = residue.atoms_dict['H']
            if atom_amide not in restraintful_set:
                restraintless_set.add(atom_amide)
    return list(restraintless_set), list(restraintful_set)


def make_proportions_plot(proteins_dict, num_bins, cs_min, cs_max):

    bin_edges, bin_midpoints = make_binning_info(num_bins, cs_min, cs_max)
    proportions_list = [[0, 0] for i in bin_midpoints]
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            protein = proteins_dict[pdb_id][bmrb_id]
            restraintless_atoms_list, restraintful_atoms_list = (
                sort_atoms_by_restraint(protein)
            )
            restraintless_binned_list = bin_atoms(
                restraintless_atoms_list, bin_edges, bin_midpoints
            )
            restraintful_binned_list = bin_atoms(
                restraintful_atoms_list, bin_edges, bin_midpoints
            )
            for i, num_restraintless in enumerate(restraintless_binned_list):
                proportions_list[i][1] += num_restraintless
                proportions_list[i][0] += restraintful_binned_list[i]

    binned_totals = [i[0] for i in proportions_list]
    proportions_list_new = [i[0] / (i[0] + i[1]) for i in proportions_list]
    proportions_list = proportions_list_new
    
    
    data_dict = {
        'proportion': [],
        'z-score': []
    }
    
    for i, prop in enumerate(proportions_list):
        data_dict['proportion'].append(prop)
        data_dict['z-score'].append(bin_midpoints[i])
    print(data_dict)
    df = pd.DataFrame(data_dict)
    
    fig = px.scatter(
        df, x='z-score', y='proportion',
        labels={
            'proportion': 'Proportion',
            'z-score': "Z(δ)"
        },
        title= "Proportion of Amide Hydrogens with NOE Restraint to Aromatic Ring Hydrogen"
    )
    fig.show(renderer="firefox")
    '''
    fig, ax = plt.subplots()
    ax.set_title("Amide-Aromatic NOEs")
    ax.set_xlabel("Secondary Shift (sigma)")
    ax.set_ylabel("Proportion of Amides w/ NOE to Aromatic Proton")
    ax.scatter(bin_midpoints, proportions_list)
    ax.set_ylim((0, 0.85))

    for i, bm in enumerate(bin_midpoints):
        ax.annotate(str(binned_totals[i]), xy=(bm, proportions_list[i] + 0.04), fontweight='bold')
    plt.show()
    '''

def make_res_prop_plot(
    proteins_dict, pairs_dict_entries, num_bins, cs_min, cs_max
):
    #first need to make restraintful_dict
    restraintful_dict = {
        'HIS': [],
        'TRP': [],
        'PHE': [],
        'TYR': [],
        'Multiple': []
    }

    restraintful_set = set()
    restraintless_set = set()
    num_double = 0
    num_triple = 0
    num_quad = 0
    num_quin = 0
    num_pairs = 0
    for pdb_id in pairs_dict_entries:
        for bmrb_id in pairs_dict_entries[pdb_id]:
            pairs_dict = pairs_dict_entries[pdb_id][bmrb_id]
            protein = proteins_dict[pdb_id][bmrb_id]
            for atom_amide in pairs_dict:
                restraintful_set.add(atom_amide)
                if len(pairs_dict[atom_amide]) == 1:
                    for res_index in pairs_dict[atom_amide]:
                        num_pairs += 1
                        atoms_aroma = pairs_dict[atom_amide][res_index]
                        res_label = atoms_aroma[0][0].res_label
                        restraintful_dict[res_label].append(atom_amide)
                else:
                    restraintful_dict['Multiple'].append(atom_amide)
                    for res_index in pairs_dict[atom_amide]:
                        num_pairs += 1
                        atoms_aroma = pairs_dict[atom_amide][res_index]
                        res_label = atoms_aroma[0][0].res_label
        
            for res_index in protein.residues_dict:
                residue = protein.residues_dict[res_index]
                if 'H' in residue.atoms_dict:
                    atom_amide = residue.atoms_dict['H']
                    if atom_amide not in restraintful_set:
                        restraintless_set.add(atom_amide)


    bin_edges, bin_midpoints = make_binning_info(num_bins, cs_min, cs_max)
    totals_dict = {
        'HIS': [0 for i in bin_midpoints],
        'TRP': [0 for i in bin_midpoints],
        'PHE': [0 for i in bin_midpoints],
        'TYR': [0 for i in bin_midpoints],
        'Multiple': [0 for i in bin_midpoints]
    }
    restraintless_binned_list = bin_atoms(
        list(restraintless_set), bin_edges, bin_midpoints
    )
    for res_label in restraintful_dict:
        atoms_list = restraintful_dict[res_label]
        totals_dict[res_label] = bin_atoms(atoms_list, bin_edges, bin_midpoints)

    binned_totals = [0 for i in bin_midpoints]
    for res_label in totals_dict:
        totals_list = totals_dict[res_label]
        for i, num in enumerate(totals_list):
            binned_totals[i] += num
    for i, num in enumerate(restraintless_binned_list):
        binned_totals[i] += num
    proportions_dict = {}
    for res_label in totals_dict:
        proportions_dict[res_label] = []
        totals_list = totals_dict[res_label]
        for i, num in enumerate(totals_list):
            proportions_dict[res_label].append(num / binned_totals[i])

    data_dict = {
        'proportion': [],
        'z-score': [],
        'res_label': []
    }
    
    for res_label in proportions_dict:
        proportions_list = proportions_dict[res_label]
        for i, prop in enumerate(proportions_list):
            data_dict['proportion'].append(prop)
            data_dict['z-score'].append(bin_midpoints[i])
            data_dict['res_label'].append(res_label)
    df = pd.DataFrame(data_dict)
    fig = px.scatter(
        df, x='z-score', y='proportion', color='res_label',
        labels={
            'proportion': 'Proportion',
            'z-score': "Z(δ)"#<sub>s</sub>
        },
        title= "Proportion of Amides w/ NOE to Aromatic Ring by Ring Type"
    )
    fig.show(renderer="firefox")














def make_res_bar_plot(pairs_dict_entries, num_bins, cs_min, cs_max):

    bin_edges, bin_midpoints = make_binning_info(num_bins, cs_min, cs_max)
    proportions_list = [[0, 0] for i in bin_midpoints]
    proportions_dict = {
        "HIS": [0 for i in bin_midpoints],
        "TRP": [0 for i in bin_midpoints],
        "PHE": [0 for i in bin_midpoints],
        "TYR": [0 for i in bin_midpoints],
    }
    for pdb_id in pairs_dict_entries:
        for bmrb_id in pairs_dict_entries[pdb_id]:
            pairs_dict = pairs_dict_entries[pdb_id][bmrb_id]
            atoms_dict = {
                "HIS": [],
                "TRP": [],
                "PHE": [],
                "TYR": []
            }
            for atom_amide in pairs_dict:
                for res_index in pairs_dict[atom_amide]:
                    atoms_aroma = pairs_dict[atom_amide][res_index]
                    if len(atoms_aroma) > 0:
                        res_label = atoms_aroma[0][0].res_label
                        atoms_dict[res_label].append(atom_amide)

            for res_label in atoms_dict:
                atoms_list = atoms_dict[res_label]
                binned_list = bin_atoms(
                    atoms_list, bin_edges, bin_midpoints
                )
                for i, num_atoms in enumerate(binned_list):
                    proportions_dict[res_label][i] += num_atoms
                    #proportions_dict[res_label][i][0] += restraintful_binned_list[i]
    
    #binned_totals = [i[0] for i in proportions_list] not bothering
    totals = [0 for i in bin_midpoints]
    for res_label in proportions_dict:
        proportions_list = proportions_dict[res_label]
        for i, num_atoms in enumerate(proportions_list):
            totals[i] += num_atoms
    
    for res_label in proportions_dict:
        proportions_list = proportions_dict[res_label]
        for i, num_atoms in enumerate(proportions_list):
            proportions_list[i] = num_atoms / totals[i]
        proportions_dict[res_label] = proportions_list
    data_dict = {
        'proportion': [],
        'z-score': [],
        'res_label': []
    }
    
    for res_label in proportions_dict:
        proportions_list = proportions_dict[res_label]
        for i, prop in enumerate(proportions_list):
            data_dict['proportion'].append(prop)
            data_dict['z-score'].append(bin_midpoints[i])
            data_dict['res_label'].append(res_label)
    df = pd.DataFrame(data_dict)
    fig = px.bar(
        df, x='z-score', y='proportion', color='res_label',
        labels={
            'proportion': 'Proportion',
            'z-score': "Z(δ)"#<sub>s</sub>
        },
        title= "Relative Proportions of Restrained Amide-Aromatic Pairs (excluding Ambiguous Restraints)"
    )
    fig.show(renderer="firefox")