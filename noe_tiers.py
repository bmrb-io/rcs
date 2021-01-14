from atoms import *

def get_confidence_tier(atoms_aroma):
    labels_aroma = [atom_aroma.atom_label for atom_aroma in atoms_aroma]
    if atoms_aroma[0].res_label == 'TRP':
        confidence_tier = get_trp_tier(labels_aroma)
    elif atoms_aroma[0].res_label == 'HIS':
        confidence_tier = get_his_tier(labels_aroma)
    elif atoms_aroma[0].res_label == 'PHE':
        confidence_tier = get_phe_tier(labels_aroma, atoms_aroma)
    elif atoms_aroma[0].res_label == 'TYR':
        confidence_tier = get_tyr_tier(labels_aroma, atoms_aroma)
    return confidence_tier, atoms_aroma[0].res_label


def get_his_tier(labels_aroma):
    high_combos = [
        ['HE1', 'HD2'], ['HE2', 'HD1']
    ]
    low_combos = [
        ['HD1', 'HE1'], ['HD2', 'HE2']
    ]
    for high_combo in high_combos:
        if high_combo[0] in labels_aroma and high_combo[1] in labels_aroma:
            return 'high'
    for low_combo in low_combos:
        if low_combo[0] in labels_aroma and low_combo[1] in labels_aroma:
            return 'low'
    return 'none'

def get_trp_tier(labels_aroma):
    high_combos = [
        ['HE3', 'HZ2'], ['HE3', 'HE1'], ['HE3', 'HD1'],
        ['HZ3', 'HZ2'], ['HZ3', 'HE1'], ['HZ3', 'HD1'],
        ['HZ2', 'HD1'], ['HZ2', 'HE1'],
        ['HH2', 'HD1'], ['HH2', 'HE1']
    ]
    low_combos = [
        ['HE1', 'HD1'], ['HE3', 'HZ3'], ['HE3', 'HH2'],
        ['HH2', 'HZ3'], ['HH2', 'HZ2']
    ]
    for high_combo in high_combos:
        if high_combo[0] in labels_aroma and high_combo[1] in labels_aroma:
            return 'high'
    for low_combo in low_combos:
        if low_combo[0] in labels_aroma and low_combo[1] in labels_aroma:
            return 'low'
    return 'none'

def get_tyr_tier(labels_aroma, atoms_aroma):
    if 'HD1' in labels_aroma and 'HD2' in labels_aroma:
        for atom in atoms_aroma:
            if atom.atom_label == 'HD1':
                cd1 = atom.cs_sigma
            elif atom.atom_label == 'HD2':
                cd2 = atom.cs_sigma
        if cd1 != cd2:
            return 'high'
        else:
            if 'QE' in labels_aroma and 'QD' in labels_aroma:
                return 'intermediate' #maybe should be high!
            else:
                return 'low'
    if 'HE1' in labels_aroma and 'HE2' in labels_aroma:
        for atom in atoms_aroma:
            if atom.atom_label == 'HE1':
                ce1 = atom.cs_sigma
            elif atom.atom_label == 'HE2':
                ce2 = atom.cs_sigma
        if ce1 != ce2:
            return 'high'
        else:
            if 'QE' in labels_aroma and 'QD' in labels_aroma:
                return 'intermediate' #maybe should be high!
            else:
                return 'low'
    if 'QE' in labels_aroma and 'QD' in labels_aroma:
        return 'intermediate'
    return 'none'

def get_phe_tier(labels_aroma, atoms_aroma):
    intermediate_combos = [
        ['QE', 'QD'], ['QE', 'HZ'], ['QD', 'HZ']
    ]
    if 'HD1' in labels_aroma and 'HD2' in labels_aroma:
        for atom in atoms_aroma:
            if atom.atom_label == 'HD1':
                cd1 = atom.cs_sigma
            elif atom.atom_label == 'HD2':
                cd2 = atom.cs_sigma
        if cd1 != cd2:
            return 'high'
        elif 'HZ' in labels_aroma:
            return 'high'
        else:
            for combo in intermediate_combos:
                if combo[0] in labels_aroma and combo[1] in labels_aroma:
                    return 'intermediate'
            return 'low'
    if 'HE1' in labels_aroma and 'HE2' in labels_aroma:
        for atom in atoms_aroma:
            if atom.atom_label == 'HE1':
                ce1 = atom.cs_sigma
            elif atom.atom_label == 'HE2':
                ce2 = atom.cs_sigma
        if ce1 != ce2:
            return 'high'
        elif 'HZ' in labels_aroma:
            return 'high'
        else:
            for combo in intermediate_combos:
                if combo[0] in labels_aroma and combo[1] in labels_aroma:
                    return 'intermediate'
            return 'low'
    for combo in intermediate_combos:
        if combo[0] in labels_aroma and combo[1] in labels_aroma:
            return 'intermediate'

    return 'none'
            


