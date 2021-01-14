
def sort_restraints_by_amide_and_ring(restraints_dict):
    sorted_dict = {}
    for restraint_id in restraints_dict:
        if len(restraints_dict[restraint_id].keys()) == 1:
            member_id = list(restraints_dict[restraint_id].keys())[0]
            restraint = restraints_dict[restraint_id][member_id]
            atom_amide = restraint.atom_amide
            atom_aroma = restraint.atom_aroma
            if atom_amide not in sorted_dict:
                sorted_dict[atom_amide] = {}
            if atom_aroma.res_index not in sorted_dict[atom_amide]:
                sorted_dict[atom_amide][atom_aroma.res_index] = []
            sorted_dict[atom_amide][atom_aroma.res_index].append(atom_aroma)
    return sorted_dict    
'''
def check_label_list(labels_list):
    d_bool = False
    e_bool = False
    if ("HD1" in labels_list and "HD2" in labels_list): #or "QD" in labels_list:
        d_bool = True
    if ("HE1" in labels_list and "HE2" in labels_list): #or "QE" in labels_list:
        e_bool = True
    return d_bool, e_bool

def check_atom_diffs(residue, d_or_e):
    #d_or_e = 'D' or 'E'
    atoms_dict = residue.atoms_dict
    atom_1 = atoms_dict['H' + d_or_e + '1']
    atom_2 = atoms_dict['H' + d_or_e + '2']
    cs_1 = atom_1.cs_sigma
    cs_2 = atom_2.cs_sigma
    if cs_1 is not None and cs_2 is not None:
        if cs_1 != cs_2:
            return True
    return False

def check_criteria(residue, labels_list):
    crit_bool = False
    d_bool, e_bool = check_label_list(labels_list)
    
    if d_bool or e_bool:
        if d_bool:
            d_or_e = 'D'
            if check_atom_diffs(residue, d_or_e):
                crit_bool = True
        if e_bool:
            d_or_e = 'E'
            if check_atom_diffs(residue, d_or_e):
                crit_bool = True
    return crit_bool, d_bool or e_bool

def find_criteria_satisfiers(protein):
    num_pairs = 0 #number of amide-ring pairs in sorted_dict
    num_satisfying_h = 0 #number with requisite two restraints
    num_w_no_shifts = 0 #number removed simply because ring proton shifts not reported
    satisfiers_list = []
    sorted_dict = sort_restraints_by_amide_and_ring(protein.restraints_dict)
    for atom_amide in sorted_dict:
        num_pairs += len(sorted_dict[atom_amide])
        for res_index in sorted_dict[atom_amide]:
            residue = protein.residues_dict[res_index]
            labels_list = [
                atom.atom_label for atom in sorted_dict[atom_amide][res_index]
            ]
            crit_bool, d_or_e_bool = check_criteria(residue, labels_list)
            if d_or_e_bool:
                num_satisfying_h += 1
            if crit_bool:
                satisfiers_list.append(atom_amide)
    return satisfiers_list, num_pairs, num_satisfying_h
'''



def check_criteria(labels_list):
    labels_list_new = []
    for label in labels_list:
        if label[0] == 'H' and len(label) > 1:
            labels_list_new.append(label)
    return len(labels_list_new) > 1

def find_criteria_satisfiers(protein):
    num_pairs = 0 #number of amide-ring pairs in sorted_dict
    satisfiers_list = []
    sorted_dict = sort_restraints_by_amide_and_ring(protein.restraints_dict)
    for atom_amide in sorted_dict:
        num_pairs += len(sorted_dict[atom_amide])
        for res_index in sorted_dict[atom_amide]:
            residue = protein.residues_dict[res_index]
            labels_list = [
                atom.atom_label for atom in sorted_dict[atom_amide][res_index]
            ]
            crit_bool = check_criteria(labels_list)
            if crit_bool:
                satisfiers_list.append(atom_amide)
    return satisfiers_list, num_pairs





