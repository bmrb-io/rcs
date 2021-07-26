import numpy as np
from atoms import Atom
from proteins import Protein
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import Union, List, Tuple, Dict

def make_binning_info(
    num_bins: int, cs_min: Union[float, int], cs_max: Union[float, int]
) -> Tuple[List[float], List[float]]:
    """
    Make lists of bin midpoints and bin edges, each in ascending order.

    Keyword arguments:
    num_bins -- number of Z-score bins to organize the atoms into for the plot
    cs_min -- the minimum Z-score to be included
    cs_max -- the maximum Z-score to be included
    Returns:
    bin_edges -- non-redundant list of the edges of the Z-score bins
    bin_midpoints -- list of the midpoints of the Z-score bins
    """
    bin_edges = list(np.linspace(cs_min, cs_max, num_bins+1))
    bin_midpoints = [
        (bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(bin_edges) - 1)
    ]
    return bin_edges, bin_midpoints

def bin_atoms(
    atoms_list: List[Atom], bin_edges: List[float], bin_midpoints: List[float]
) -> List[int]:
    """
    Count number of atoms to be placed into each z-score bin.

    Keyword arguments:
    atoms_list -- list of Atom instances
    bin_edges -- non-redundant list of the edges of the Z-score bins
    bin_midpoints -- list of the midpoints of the Z-score bins
    Returns:
    binned_list -- list containing number of atoms in each bin (according to
        index)
    """
    cs_list = []
    for atom in atoms_list:
        cs_sigma = float(atom.cs_sigma)
        if abs(cs_sigma) <= bin_edges[-1]: #Right now only works if cs_min = -cs_max
            cs_list.append(cs_sigma)
    inds = np.digitize(cs_list, bin_edges)
    binned_list = [0 for i in bin_midpoints]
    for bin_ind in inds:
        bin_ind = bin_ind - 1
        binned_list[bin_ind] += 1
    return binned_list

def sort_atoms_by_restraint(protein: Protein) -> Tuple[List[Atom], List[Atom]]:
    """
    Separate Atoms in a Protein with a restraint to an aromatic ring proton 
    from those without.

    Keyword arguments:
    protein -- a Protein instance
    Returns:
    restraintless_list -- list of Atoms from Protein without a restraint
    restraintful_list -- list of Atoms from Protein with a restraint
    """
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
    restraintless_list = list(restraintless_set)
    restraintful_list = list(restraintful_set)
    return restraintless_list, restraintful_list


def make_proportions_plot(
    proteins_dict: Dict[str, Dict[str, Protein]], num_bins: int, 
    cs_min: Union[float, int], cs_max: Union[float, int], fig
):
    """
    Generate a plot of the proportion of amides from proteins_dict that have
    at least one restraint to an aromatic-ring proton for varying values of 
    the amide Z-score.

    Keyword arguments:
    proteins_dict -- dict of Protein instances organized by PDB and BMRB ID
    num_bins -- number of Z-score bins to organize the atoms into for the plot
    cs_min -- the minimum Z-score to be included
    cs_max -- the maximum Z-score to be included
    """
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

    binned_totals = [i[0] + i[1] for i in proportions_list] #total atoms in this bin
    proportions_list_new = [i[0] / (i[0] + i[1]) for i in proportions_list] #proportion with restraints
    proportions_list = proportions_list_new

    fig.add_trace(
        go.Scatter(
            mode='markers+text',
            x=bin_midpoints,
            y=proportions_list,
            marker=dict(
                color='#7F7F7F',
                size=20,
                line=dict(width=2)
            ),
            text=[str(i) for i in binned_totals], # Write total for that bin above marker
            textposition='top center',
            textfont_size=18,
            textfont_family="Arial",
            showlegend=False,
        ), row=1, col=1
    )
    return fig

def make_num_restraints_plot(
    proteins_dict: Dict[str, Dict[str, Protein]], 
    outlier_sigma: Union[int, float]
):
    num_restraints_dict = {
        'Upfield': {'HIS': {}, 'TRP': {}, 'PHE': {}, 'TYR': {}},
        'Downfield': {'HIS': {}, 'TRP': {}, 'PHE': {}, 'TYR': {}},
        'Normal': {'HIS': {}, 'TRP': {}, 'PHE': {}, 'TYR': {}},
    }
    num_restraints_dict = {
        'HIS': {
            '1': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            '2': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            #'>1': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            '>2': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
        },
        'TRP': {
            '1': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            #'>1': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            '2': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            '>2': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
        },
        'PHE': {
            '1': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            #'>1': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            '2': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            '>2': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
        },
        'TYR': {
            '1': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            #'>1': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            '2': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
            '>2': {'Upfield': 0, 'Normal': 0, 'Downfield': 0},
        },
    }
    num_redun = 0
    redun_set = set()
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            protein = proteins_dict[pdb_id][bmrb_id]
            protein.prune_undefined_pairs()
            pairs_dict = protein.pairs_dict
            for atom_amide in pairs_dict:
                if atom_amide.cs_sigma >= outlier_sigma:
                    shift_type = 'Downfield'
                elif atom_amide.cs_sigma <= -1 * outlier_sigma:
                    shift_type = 'Upfield'
                else:
                    shift_type = 'Normal'
                if len(pairs_dict[atom_amide]) == 1:
                    for res_index_aroma in pairs_dict[atom_amide]: #only one, but simpler this way
                        atoms_aroma = pairs_dict[atom_amide][res_index_aroma]
                        labels_aroma = [atom[0].atom_label for atom in atoms_aroma]
                        num_restraints = len(set(labels_aroma))
                        if num_restraints > 0:
                            if num_restraints > 2:
                                num_restraints = '>2'
                            #if num_restraints > 1:
                            #    num_restraints = '>1'
                            else:
                                num_restraints = str(num_restraints)
                            res_label_aroma = atoms_aroma[0][0].res_label
                            num_restraints_dict[res_label_aroma][num_restraints][shift_type] += 1
    totals_by_res = {
        'HIS': {},
        'TRP': {},
        'PHE': {},
        'TYR': {}
    }
    for res_label in num_restraints_dict:
        nr_res = num_restraints_dict[res_label]
        totals = {'Upfield': 0, 'Normal': 0, 'Downfield': 0}
        for nr in nr_res:
            for shift_type in nr_res[nr]:
                num_atoms = nr_res[nr][shift_type]
                totals[shift_type] += num_atoms
        totals_by_res[res_label] = totals
    nr_dict_normalized = {}
    for res_label in num_restraints_dict:
        nr_dict_normalized[res_label] = {}
        for nr in num_restraints_dict[res_label]:
            nr_dict_normalized[res_label][nr] = {}
            for shift_type in num_restraints_dict[res_label][nr]:
                num_atoms = num_restraints_dict[res_label][nr][shift_type]
                total = totals_by_res[res_label][shift_type]
                nr_dict_normalized[res_label][nr][shift_type] = num_atoms

    fig = make_subplots(
        rows=2, cols=2, 
        subplot_titles=("<b>HIS<b>", "<b>TRP<b>", "<b>PHE<b>", "<b>TYR<b>"),
        shared_yaxes=True, vertical_spacing=0.10
    )
    for i in range(4):
        fig.layout.annotations[i].update(font=dict(family="Arial",size=18))
    i = 0
    row_nums = [1, 1, 2, 2]
    col_nums = [1, 2, 1, 2]
    colors = ['rgb(102, 197, 204)', 'rgb(248, 156, 116)', 'rgb(220, 176, 242)']

    fig.update_yaxes(row=1, col=1, title_text='Proportion of Pairs', type='log')
    fig.update_yaxes(row=1, col=2, type='log')
    fig.update_yaxes(row=2, col=1, title_text='Proportion of Pairs', type='log')
    fig.update_yaxes(row=2, col=2, type='log')

    for res_label in nr_dict_normalized:

        row_num = row_nums[i]
        col_num = col_nums[i]
        nr_dict = nr_dict_normalized[res_label]
        j = 0
        for num_restraints in nr_dict:
            totals = list(num_restraints_dict[res_label][num_restraints].values())
            restraints_by_shift = nr_dict[num_restraints]
            if i == 1:
                fig.add_trace(
                    go.Bar(
                        x=list(restraints_by_shift.keys()), 
                        y=list(restraints_by_shift.values()),
                        name=num_restraints,
                        marker_color=colors[j],
                        showlegend=True,
                        text=totals,
                        textposition='inside',
                        textfont_size=18,
                        textfont_family="Arial",
                    ),
                    row=row_num, col=col_num
                )
            else:
                fig.add_trace(
                    go.Bar(
                        x=list(restraints_by_shift.keys()), 
                        y=list(restraints_by_shift.values()),
                        name=num_restraints,
                        marker_color=colors[j],
                        showlegend=False,
                        text=totals,
                        textposition='inside',
                        textfont_size=18,
                        textfont_family="Arial",
                    ),
                    row=row_num, col=col_num
                )

            j+=1
        i+=1

    fig.update_layout(
        #title=(
        #    'Restrained Amide-Aromatic Pairs by'
        #    + '<br>'
        #    + 'Number of Restraints'
        #),
        title_x=0.5,
        title_y=0.97,
        font=dict(family="Arial",size=18),
        legend_title="Number of" + "<br>" + "Restraints",
        legend_x = 0.01,
        legend_y = 0.99,
        autosize=False,
        width=1100,
        height=1100,
        barmode='group'#'stack'
    )
    #fig.update_traces(textangle=0)

    #fig.update_yaxes(row=1, col=1, range=[0, 1.1])
    #fig.update_yaxes(row=1, col=2, range=[0, 1.1])
    #fig.update_yaxes(row=2, col=1, range=[0, 1.1])
    #fig.update_yaxes(row=2, col=2, range=[0, 1.1])

    #fig.show(renderer="firefox")
    fig.write_image("../images/noes_by_num.pdf")
    



def make_res_prop_plot(
    proteins_dict: Dict[str, Dict[str, Protein]], num_bins: int, 
    cs_min: Union[float, int], cs_max: Union[float, int], fig
):
    """
    Generate a plot of the proportion of amides from proteins_dict that have
    at least one restraint to an aromatic-ring proton for varying values of 
    the amide Z-score, and separate by res_label of the ring.

    Keyword arguments:
    proteins_dict -- dict of Protein instances organized by PDB and BMRB ID
    num_bins -- number of Z-score bins to organize the atoms into for the plot
    cs_min -- the minimum Z-score to be included
    cs_max -- the maximum Z-score to be included
    """
    restraintful_dict = {
        'HIS': [],
        'TRP': [],
        'PHE': [],
        'TYR': [],
        'Multiple': []
    }

    restraintful_set = set()
    restraintless_set = set()
    for pdb_id in proteins_dict:
        for bmrb_id in proteins_dict[pdb_id]:
            protein = proteins_dict[pdb_id][bmrb_id]
            pairs_dict = protein.pairs_dict
            for atom_amide in pairs_dict: # if it's in pairs_dict, it must have a restraint
                restraintful_set.add(atom_amide)
                if len(pairs_dict[atom_amide]) == 1:
                    for res_index in pairs_dict[atom_amide]: 
                        atoms_aroma = pairs_dict[atom_amide][res_index]
                        res_label = atoms_aroma[0][0].res_label # The res_label to assign the amide
                        restraintful_dict[res_label].append(atom_amide)
                else: # The amide has restraints to ring protons from more than one residue
                    restraintful_dict['Multiple'].append(atom_amide)
                    for res_index in pairs_dict[atom_amide]:
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
                showlegend=True,
            ), row=2,col=1
        )
    return fig

def make_all_plots(proteins_dict: Dict[str, Dict[str, Protein]], num_bins: int, 
    cs_min: Union[float, int], cs_max: Union[float, int]
):
    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.1, 
        subplot_titles=('<b>A<b>', '<b>B<b>')
    )
    fig = make_proportions_plot(proteins_dict, num_bins, cs_min, cs_max, fig)
    fig = make_res_prop_plot(proteins_dict, num_bins, cs_min, cs_max, fig)
    fig.update_yaxes(title_text='Proportion', row=1, col=1, range=[0, 0.8])
    fig.update_yaxes(title_text='Proportion', row=2, col=1, range=[0, 0.3])
    fig.update_layout(
        xaxis_title='Z-score',
        xaxis_anchor = 'y2',
        font=dict(family="Arial",size=18),
        legend=dict(
            y=0.28,
            x=0.86,
        ),
        width=1100,
        height=1100
    )
    fig.update_xaxes(title_standoff=35) 
    fig.layout.annotations[0].update(x=0.02, y=0.97, font=dict(family="Arial",size=18)) # 0.07
    fig.layout.annotations[1].update(x=0.02, y=0.42, font=dict(family="Arial",size=18)) # 0.11
    #fig.show(renderer="firefox")
    fig.write_image("../images/combo_plot.pdf")

