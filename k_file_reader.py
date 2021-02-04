from atoms import Atom
from residues import Residue
from proteins import Protein

def make_protein_from_file(filename):
    '''
    Build protein complete with amide and aromatic residues and chemical 
    shifts for amide and aromatic ring protons.

    Keyword arguments:
    filename -- path to file containing chemical shift data
    Returns:
    proteins_dict -- dict of Protein objects organize by PDB ID and BMRB ID;
        Proteins contain residues_dict but NoneType restraints_dict
    exceptions_map_entries -- dict of exceptions raised; organized by PDB ID
        and BMRB ID
    '''
    with open(filename) as infile:
        lines = infile.readlines()
        if len(lines) == 1: # A one-line exception was written to the file
            exc_str = lines[0]
            split_index = exc_str.find(' ')
            exc_str = exc_str[split_index+1:] # Getting rid of superfluous info
            return exc_str
        line0 = lines[0].split(',')
        pdb_id = line0[0]
        bmrb_id = line0[1]
        protein = Protein(pdb_id, bmrb_id)
        if line0[6] == '1' and line0[7] == '1': # There's only 1 entity/assembly
            for line in lines:
                line = line.split(',')
                protein.residues_dict, protein.exceptions_map_residues = (
                    add_residues(
                        line, protein.residues_dict, 
                        protein.exceptions_map_residues
                    )
                )
                protein.pair_geometries = add_pair_geometries(
                    line, protein.pair_geometries)
        else: #Return an exception
            return "Too many entities/assemblies"

    return protein

def add_pair_geometries(line, pair_geometries):
    """
    Add distances between the amide and the aromatic ring protons listed on 
    a line of the k-file to the pair_geometries.

    Keyword arguments:
    line -- list of the data values on a line of the k-file
    pair_geometries -- dict organized by res_index of amide and aromatic, then
        atom_label of aromatic ring proton; values are distance between atoms
    Returns:
    pair_geometries -- with the distance values for 5 nearest aromatic rings
        and their protons to the amide proton on the line
    """
    atom_labels_in_file = {
        "PHE": ['HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
        "TYR": ['HD1', 'HD2', 'HE1', 'HE2', 'HH'],
        "TRP": ['HE3', 'HZ2', 'HZ3', 'HH2', 'HE1'],
        "HIS": ['HD2', 'HE1', 'HE2', 'xx', 'yy']
    } #order of the atoms in the file
    num = 52
    res_index_amide = line[2]
    amide_geometries = {} # pair geometries for the amide on this line
    for i in range(5):
        start_index = 8 + num * i #to begin with the first ring proton 
        res_index_aroma = line[start_index]
        if res_index_aroma == ".":
            break
        amide_geometries[res_index_aroma] = {}
        res_label_aroma = line[start_index+1]
        for j in range(5):
            dist_index = start_index + 34 + j * 4 #move to next proton's distance value
            dist = float(line[dist_index])
            atom_label_aroma = atom_labels_in_file[res_label_aroma][j]
            amide_geometries[res_index_aroma][atom_label_aroma] = dist
    pair_geometries[res_index_amide] = amide_geometries
    return pair_geometries

def add_residues(line, residues_dict, exceptions_map):
    """
    Add residues from a line of the data file to residues_dict if not already
    in it; otherwise update atoms_dict of the residues.

    Keyword arguments:
    line -- a line from the data file
    residues_dict -- dict containing Residues for a PDB entry organized by 
        residue index
    exceptions_map -- dict containing exceptions generated when trying to make
        residues; organized by residue index
    """
    res_amide, res_aroma_list = make_residues(line)
    if res_amide.res_index in residues_dict: # The Residue was already made
        residues_dict[res_amide.res_index].atoms_dict.update(
            res_amide.atoms_dict
        )
    else: # This is the first appearance of the Residue in the k-file
        residues_dict[res_amide.res_index] = res_amide
    for res_aroma in res_aroma_list: # A line of the k-file includes 5 nearest rings
        if res_aroma.res_index in residues_dict:
            residues_dict[res_aroma.res_index].atoms_dict.update(
                res_aroma.atoms_dict
            )
        else:
            residues_dict[res_aroma.res_index] = res_aroma
    return residues_dict, exceptions_map

def make_residues(line):
    """
    Make one amide and multiple aromatic residues from line of ring data file.

    Keyword arguments:
    line -- a line from the data file
    Returns:
    res_amide -- Residue object containing the amide proton
    res_aroma_list -- list of Residue objects containing aromatic ring protons
    """
    res_amide = make_res_amide(line)
    res_aroma_list = []
    for i in range(5): # A line of the k-file includes 5 nearest rings
        res_aroma = make_res_aroma(line, i)
        if isinstance(res_aroma, Residue):
            res_aroma_list.append(make_res_aroma(line, i))
    return res_amide, res_aroma_list

def make_res_amide(line):
    """
    Make a residue containing the amide atom from a line of ring data file.
    
    Keyword arguments:
    line -- a line from the data file
    i -- index of the ring in the line
    Returns:
    res_amide -- Residue object
    """
    res_index = line[2]
    res_label = line[3]
    cs_sigma = float(line[5])
    atom = Atom(res_index, res_label, 'H', cs_sigma)
    atoms_dict = {atom.atom_label: atom}
    res_amide = Residue(res_index, res_label, atoms_dict)
    return res_amide

def make_res_aroma(line, i):
    """
    Make an aromatic residue from a line of Kumaran's ring data file.

    Keyword arguments:
    line -- a line from the data file
    i -- index of the ring in the line
    Returns:
    'Not enough rings' -- if data for this ring is empty
    res -- Residue object
    """
    num = 52 #the periodicity of aromatics in file
    start_index = 8 + num * i
    end_index = 8 + num + num * i
    ring_data = line[start_index:end_index]
    res_index = ring_data[0]
    res_label = ring_data[1]
    if res_index == '.': # A . gets put in for the last rings in k-file if <5 were found
        return "Not enough rings"
    atoms_list = make_atoms_aroma(ring_data)
    atoms_dict = {}
    for atom in atoms_list:
        atoms_dict[atom.atom_label] = atom
    res = Residue(res_index, res_label, atoms_dict)
    return res


def make_atoms_aroma(ring_data):
    """
    Make a list of aromatic ring protons from Kumaran's ring data file.

    Keyword arguments:
    ring_data -- a slice of a line from the data file corresponding to the 
        particular ring
    Returns:
    atoms_list -- list of Atom objects for all aromatic ring protons
    """
    atoms_file_dict = {
        "PHE": {20: 'HD1', 22: 'HD2', 24: 'HE1', 26: 'HE2', 28: 'HZ'},
        "TYR": {20: 'HD1', 22: 'HD2', 24: 'HE1', 26: 'HE2', 28: 'HH'},
        "TRP": {20: 'HE3', 22: 'HZ2', 24: 'HZ3', 26: 'HH2', 28: 'HE1'}, 
        "HIS": {20: 'HD2', 22: 'HE1', 24: 'HE2', 26: 'xx', 28: 'yy'},
    } # xx and yy are pseudoatoms to maintain the formatting in the k-file
    res_index = ring_data[0]
    res_label = ring_data[1]
    atoms_list = []
    for i in range(5): # Go through all 5 rings on the line. 
        cs_sigma_index = 20 + (2 * i)
        cs_sigma = ring_data[cs_sigma_index]
        if cs_sigma == '.':
            cs_sigma = None
        else:
            cs_sigma = float(cs_sigma)
        atom_label = atoms_file_dict[res_label][cs_sigma_index]
        atom = Atom(res_index, res_label, atom_label, cs_sigma)
        atoms_list.append(atom)
    return atoms_list

