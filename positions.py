import os
from pdbecif.mmcif_io import CifFileReader

def get_positions(pdb_id):

    file_path = './cifs/{}.cif'.format(pdb_id)
    if not check_file_exists(file_path):
        get_file(pdb_id, file_path)
    atom_site_obj, num_atoms = make_site_obj(file_path, pdb_id)
    atoms_dict = make_atoms_dict(atom_site_obj, num_atoms) 
    return atoms_dict

def check_file_exists(file_path):
    
    file_check = os.path.isfile(file_path)

    return file_check

def get_file(pdb_id, file_path):

    cmd = 'wget https://files.rcsb.org/download/{}.cif -O' + file_path
    cmd = cmd.format(pdb_id)
    os.system(cmd)

def make_site_obj(file_path, pdb_id):
    
    cfr = CifFileReader()
    cif_obj = cfr.read(file_path, output='cif_wrapper')
    cif_data = cif_obj[pdb_id]
    atom_site_obj = cif_data._atom_site

    num_models = int(atom_site_obj.pdbx_PDB_model_num[-1]) 
    num_atoms = int(
        len(atom_site_obj.pdbx_PDB_model_num) / num_models
    )

    return atom_site_obj, num_atoms

def make_atoms_dict(atom_site_obj, num_atoms):

    atoms_dict = {}

    for i in range(num_atoms):
        amino_index = atom_site_obj.label_seq_id[i]
        amino_label = atom_site_obj.label_comp_id[i]
        amino_key = (amino_index, amino_label)
        atom_label = atom_site_obj.label_atom_id[i]
        atom_position = [
            atom_site_obj.Cartn_x[i],
            atom_site_obj.Cartn_y[i],
            atom_site_obj.Cartn_z[i]
        ]
        if amino_key not in atoms_dict:
            atoms_dict[amino_key] = {}
        atoms_dict[amino_key][atom_label] = [atom_position, None]
    
    return atoms_dict
        

