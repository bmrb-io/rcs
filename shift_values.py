import pynmrstar

def get_shifts_list(bmrb_id):

    atoms_dict = {}
    entry = pynmrstar.Entry.from_database(bmrb_id)
    try:
        shifts_loop = entry.get_loops_by_category("Atom_chem_shift")[0]
        shifts_list = shifts_loop.get_tag(
            ['Comp_index_ID', 'Comp_ID', 'Atom_ID', 'Val', 'Val_err']
        )
        return shifts_list

    except IndexError:
        print(bmrb_id + " doesn't have any chemical shift values.")
        return None

def match_shifts(shifts_list, atoms_dict):
    
    nt = 0
    nf = 0

    framed_bool = False

    if not check_alignment(shifts_list, atoms_dict):
        frameshift_index = find_frameshift_index(shifts_list, atoms_dict)
        if frameshift_index is not None:
            shifts_list = shifts_list[frameshift_index:]
            framed_bool = True

    for shift in shifts_list:
        amino_key = (shift[0], shift[1])
        if amino_key in atoms_dict:
            nt+=1
        else:
            nf +=1

    return nt, nf, framed_bool

def check_alignment(shifts_list, atoms_dict, num_to_check = 5):

    for shift in shifts_list[:5]:
        amino_key = (shift[0], shift[1])
        if amino_key not in atoms_dict:
            return False
    
    return True

def find_frameshift_index(shifts_list, atoms_dict):
    
    bmrb_aminos_list, pdb_aminos_list = prepare_aminos_lists(
        shifts_list, atoms_dict
    )
    amino_index = match_frames(bmrb_aminos_list, pdb_aminos_list)
    if amino_index is not None:
        for i, shift in enumerate(shifts_list):
            if shift[0] == str(amino_index + 1):
                frameshift_index = i
                return frameshift_index
    else:
        return None



def prepare_aminos_lists(shifts_list, atoms_dict):

    pdb_aminos_list = list(atoms_dict.keys())
    pdb_aminos_list_temp = []
    for amino in pdb_aminos_list:
        if amino[0] == '.':
            pdb_aminos_list.remove(amino)
        else:
            pdb_aminos_list_temp.append((int(amino[0]), amino[1]))
            
    pdb_aminos_list = pdb_aminos_list_temp
    pdb_aminos_list.sort()
    bmrb_aminos_list = list(set([(int(shift[0]), shift[1]) for shift in shifts_list]))
    bmrb_aminos_list.sort()

    return bmrb_aminos_list, pdb_aminos_list

def match_frames(bmrb_aminos_list, pdb_aminos_list):

    pdb_amino_index, pdb_amino_label = pdb_aminos_list[0]
    match_count = 0
    for i, bmrb_amino in enumerate(bmrb_aminos_list):
        bmrb_amino_label = bmrb_amino
        for k in range(5):
            if bmrb_aminos_list[i+k][1] == pdb_aminos_list[k][1]:
                match_count += 1
                if match_count == 5:
                    return i
            else:
                break

    return None
