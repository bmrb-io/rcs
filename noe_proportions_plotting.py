import numpy as np
import plotly.graph_objects as go

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

    binned_totals = [i[0] + i[1] for i in proportions_list]
    proportions_list_new = [i[0] / (i[0] + i[1]) for i in proportions_list]
    proportions_list = proportions_list_new


    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            mode='markers+text',
            x=bin_midpoints,
            y=proportions_list,
            marker=dict(
                color='LightSkyBlue',
                size=20,
                line=dict(
                    color='MediumPurple',
                    width=2
                )
            ),
            text=[str(i) for i in binned_totals],
            textposition='top center',
            textfont_size=16,
            textfont_family="Courier New, monospace",
            showlegend=False
        )
    )
    fig.update_layout(
        title=(
            'Proportion of Amide Hydrogens with NOE Restraint'
            + '<br>'
            + 'to Aromatic Ring Hydrogen'
        ),
        title_x=0.5,
        xaxis_title='Z(δ)',
        yaxis_title='Proportion',
        font=dict(family="Courier New, monospace",size=16),
        yaxis_range=[0,0.78]
    )
    fig.show(renderer="firefox")


def make_res_prop_plot(
    proteins_dict, num_bins, cs_min, cs_max
):
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
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            protein = proteins_dict[pdb_id][bmrb_id]
            pairs_dict = protein.pairs_dict
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
    fig = go.Figure()
    for res_label in proportions_dict:
        proportions_list = proportions_dict[res_label]
        totals_list = totals_dict[res_label]
        fig.add_trace(
            go.Scatter(
                mode='markers',
                name=res_label,
                x=bin_midpoints,
                y=proportions_list,
                marker=dict(
                    size=20,
                    line=dict(width=2)
                ),
                showlegend=True
            )
        )
    fig.update_layout(
        title=(
            'Proportion of Amide Hydrogens with NOE Restraint'
            + '<br>'
            + 'to Aromatic Ring Hydrogen'
        ),
        title_x=0.5,
        xaxis_title='Z(δ)',
        yaxis_title='Proportion',
        font=dict(family="Courier New, monospace",size=16),
        yaxis_range=[0,0.35]
    )
    fig.show(renderer="firefox")