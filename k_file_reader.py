from atoms import Atom
from residues import Residue
from proteins import Protein

def make_proteins_from_file(filename):
    exceptions_map_entries = {}
    proteins_dict = {}
    with open(filename) as infile:
        lines = infile.readlines()[157:] #should be 157:
        for line in lines:
            line = line.split(',')
            if line[6] == '1' and line[7] == '1':
                pdb_id = line[0]
                bmrb_id = line[1]
                if pdb_id not in proteins_dict:
                    proteins_dict[pdb_id] = {}
                if bmrb_id not in proteins_dict[pdb_id]:
                    proteins_dict[pdb_id][bmrb_id] = Protein(pdb_id, bmrb_id)
                protein = proteins_dict[pdb_id][bmrb_id]
                protein.residues_dict, protein.exceptions_map_residues = (
                    add_residues(
                        line, protein.residues_dict, 
                        protein.exceptions_map_residues
                    )
                )
                proteins_dict[pdb_id][bmrb_id] = protein
            else:
                pdb_id = line[0]
                bmrb_id = line[1]
                if pdb_id not in exceptions_map_entries:
                    exceptions_map_entries[pdb_id] = {}
                exceptions_map_entries[pdb_id][bmrb_id] = (
                    "Too many entities/assemblies"
                )
    return proteins_dict, exceptions_map_entries

def add_residues(line, residues_dict, exceptions_map):
    res_amide, res_aroma_list = make_residues(line)
    if res_amide.res_index in residues_dict:
        residues_dict[res_amide.res_index].atoms_dict.update(
            res_amide.atoms_dict
        )
    else:
        residues_dict[res_amide.res_index] = res_amide
    for res_aroma in res_aroma_list:
        if res_aroma.res_index in residues_dict:
            residues_dict[res_aroma.res_index].atoms_dict.update(
                res_aroma.atoms_dict
            )
        else:
            residues_dict[res_aroma.res_index] = res_aroma
    return residues_dict, exceptions_map

def make_residues(line):
    res_amide = make_res_amide(line)
    res_aroma_list = []
    for i in range(5):
        res_aroma = make_res_aroma(line, i)
        if isinstance(res_aroma, Residue):
            res_aroma_list.append(make_res_aroma(line, i))
    return res_amide, res_aroma_list

def make_res_amide(line):
    res_index = line[2]
    res_label = line[3]
    cs_sigma = float(line[5])
    atom = Atom(res_index, res_label, 'H', cs_sigma)
    atoms_dict = {atom.atom_label: atom}
    res_amide = Residue(res_index, res_label, atoms_dict)
    return res_amide

def make_res_aroma(line, i):
    num = 30 #the periodicity of aromatics in file
    start_index = 8 + num * i
    end_index = 8 + num + num * i
    ring_data = line[start_index:end_index]
    res_index = ring_data[0]
    res_label = ring_data[1]
    if res_index == '.':
        return "Not enough rings"
    atoms_list = make_atoms_aroma(ring_data)
    atoms_dict = {}
    for atom in atoms_list:
        atoms_dict[atom.atom_label] = atom
    res = Residue(res_index, res_label, atoms_dict)
    return res


def make_atoms_aroma(ring_data):
    atoms_file_dict = {
        "PHE": {20: 'HD1', 22: 'HD2', 24: 'HE1', 26: 'HE2', 28: 'HZ'},
        "TYR": {20: 'HD1', 22: 'HD2', 24: 'HE1', 26: 'HE2', 28: 'HH'},
        "TRP": {20: 'HE3', 22: 'HZ2', 24: 'HZ3', 26: 'HH2', 28: 'HE1'}, 
        "HIS": {20: 'HD2', 22: 'HE1', 24: 'HE2', 26: 'xx', 28: 'yy'}, #this is messed up, but works for now
    }
    res_index = ring_data[0]
    res_label = ring_data[1]
    atoms_list = []
    for i in range(5):
        cs_sigma_index = 20 + (2 * i)
        #ambi_index = 29 + (2 * i)
        ##ignoring ambiguity codes for now!
        cs_sigma = ring_data[cs_sigma_index]
        if cs_sigma == '.':
            cs_sigma = None
        else:
            cs_sigma = float(cs_sigma)
        atom_label = atoms_file_dict[res_label][cs_sigma_index]
        atom = Atom(res_index, res_label, atom_label, cs_sigma)
        atoms_list.append(atom)
    return atoms_list

