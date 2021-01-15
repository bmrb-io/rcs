from atoms import Atom
from residues import Residue
from restraints import Restraint
from proteins import Protein
from pdbecif.mmcif_io import CifFileReader

def get_coordinates(pdb_id):
    """
    Extract coordinate information from cif file to build protein.
    """
    file_path = f"./cifs/{pdb_id}.cif"
    if not check_file_exists(file_path)
        get_file(pdb_id, file_path)
    atom_site_obj = make_site_obj(file_path, pdb_id)


def check_file_exists(file_path):
    "Check if there is a cif file at expected location."
    file_check = os.path.isfile(file_path)    
    return file_check

def get_file(pdb_id, file_path):
    "Download cif file from RCSB and move to file_path."
    cmd = f"wget https://files.rcsb.org/download/{pdb_id}.cif -O" + file_path
    os.system(cmd)

def make_site_obj(file_path, pdb_id):
    cfr = CifFileReader()
    cif_obj = cfr.read(file_path, output='cif_wrapper')
    cif_data = cif_obj[pdb_id]
    atom_site_obj = cif_data._atom_site

def make_residues_dict(atom_site_obj):
    residues_dict = {}
    for i in range(num_atoms):
        res_index = atom_site_obj.labe_seq_id[i]
        if res_index not in residues_dict:
            residues_dict[res_index]['res_label'] = res_label
        atom_label = atom_site_obj.label_atom_id[i]