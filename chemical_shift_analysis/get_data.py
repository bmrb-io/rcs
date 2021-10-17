import pynmrstar
from mmcif.io.PdbxReader import PdbxReader
import os
import numpy
import gzip
from typing import Union, List, Optional

aromatic_atoms = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
    'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1'],
    'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', 'HD1', 'HD2', 'HE1', 'HE2', 'xx', 'yy']  # if needed un comment
}


def get_chemical_shifts(str_file: str,auth_tag: Optional[bool] =False) -> tuple:
    """
    Extract the chemical shift information from  the NMR-STAR file
    :param str_file: NMR-STAR file name with full path
    :param auth_tag: use sequence id from Auth_seq_ID, optional
    :return: tuple of amide chemical  shift as dictionary, aromatic chemical shift as dictionary, entity size
        and entity and entity assembly size
    """
    str_data = pynmrstar.Entry.from_file(str_file)
    csdata = str_data.get_loops_by_category('Atom_chem_shift')
    entity = str_data.get_tag('_Entity_assembly.Entity_ID')
    entity_size = len(set(entity))
    assembly_size = (len(entity))
    amide_cs = {}
    aromatic_cs = {}
    for cs in csdata:
        tag_list = cs.get_tag_names()

        id2 = tag_list.index('_Atom_chem_shift.Auth_asym_ID')
        if auth_tag:
            id1 = tag_list.index('_Atom_chem_shift.Auth_seq_ID')
        else:
            id1 = tag_list.index('_Atom_chem_shift.Comp_index_ID')
        id3 = tag_list.index('_Atom_chem_shift.Comp_ID')
        id4 = tag_list.index('_Atom_chem_shift.Atom_ID')
        id5 = tag_list.index('_Atom_chem_shift.Val')
        id6 = tag_list.index('_Atom_chem_shift.Ambiguity_code')
        for d in cs.data:
            if d[id4] == 'H':
                if d[id2] == '.':
                    d[id2] = 'A'  # temp fix if asym id is missing
                amide_cs[(d[id1], d[id2], d[id3], d[id4])] = (d[id5], get_z_score(d[id3], float(d[id5])))
            if d[id3] in aromatic_atoms.keys():
                if d[id2] == '.':
                    d[id2] = 'A'  # temp fix if asym id is missing
                k = (d[id1], d[id2], d[id3])
                if d[id4] in aromatic_atoms[d[id3]]:
                    if k not in aromatic_cs.keys():
                        aromatic_cs[k] = {}
                    aromatic_cs[k][d[id4]] = (d[id5], d[id6])
    return amide_cs, aromatic_cs, entity_size, assembly_size


def get_coordinates(cif_file: str, use_auth_tag: Optional[bool] =True, nmrbox: Optional[bool] =False) -> dict:
    """
    Extract coordinate information from cif file as a dictionary
    {model_id : {(seq_id,chain_id,res_id,atom_id) : array[x,y,x],...},...}
    :param cif_file: Input CIF file with full path
    :param use_auth_tag: use sequence id from author provided tag
    :param nmrbox: use the gzip CIF file from NMRBox reboxitory
    :return: dictionary
    """
    cif_data = []
    if nmrbox:
        ifh = gzip.open(cif_file,'rt')
    else:
        ifh = open(cif_file, 'r')
    pRd = PdbxReader(ifh)
    pRd.read(cif_data)
    ifh.close()
    c0 = cif_data[0]
    atom_site = c0.getObj('atom_site')
    max_models = int(atom_site.getValue('pdbx_PDB_model_num', -1))
    col_names = atom_site.getAttributeList()
    model_id = col_names.index('pdbx_PDB_model_num')
    x_id = col_names.index('Cartn_x')
    y_id = col_names.index('Cartn_y')
    z_id = col_names.index('Cartn_z')
    atom_id = col_names.index('label_atom_id')
    comp_id = col_names.index('label_comp_id')
    asym_id = col_names.index('label_asym_id')
    entity_id = col_names.index('label_entity_id')
    seq_id = col_names.index('label_seq_id')
    icode_id = col_names.index('pdbx_PDB_ins_code')
    alt_id = col_names.index('label_alt_id')
    aut_seq_id = col_names.index('auth_seq_id')
    aut_asym_id = col_names.index('auth_asym_id')
    aut_atom_id = col_names.index('auth_atom_id')
    aut_comp_id = col_names.index('auth_comp_id')
    pdb_models = {}
    atom_ids = {}
    for model in range(1, max_models + 1):
        pdb = {}
        aid = {}
        for dat in atom_site.getRowList():
            if dat[atom_id] == 'H' or dat[comp_id] in aromatic_atoms.keys():  # Only necessary coordinates for this
                # calculation
                if int(dat[model_id]) == model:
                    if use_auth_tag:
                        aid[(dat[aut_seq_id], dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
                            (dat[entity_id], dat[asym_id], dat[comp_id], dat[seq_id], dat[aut_seq_id],
                             dat[alt_id], dat[icode_id], dat[aut_asym_id])
                        pdb[(dat[aut_seq_id], dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
                            numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
                    else:
                        aid[(dat[seq_id], dat[asym_id], dat[comp_id], dat[atom_id])] = \
                            (dat[entity_id], dat[asym_id], dat[comp_id], dat[seq_id], dat[aut_seq_id],
                             dat[alt_id], dat[icode_id], dat[aut_asym_id])
                        pdb[(dat[seq_id], dat[asym_id], dat[comp_id], dat[atom_id])] = \
                            numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
        pdb_models[model] = pdb
        atom_ids[model] = aid
    return pdb_models


def get_pdb_data(pdb_id: str, auth_tag: Optional[bool] =False, nmrbox: Optional[bool]=False):
    """
    Extract coordinate data as dictionary for a given PDB ID
    :param pdb_id: PDB ID
    :param auth_tag: use sequence information from author tag
    :param nmrbox: Instead of download, use the gzip CIF file from NMRBox reboxitory
    :return: dictionary
    """
    if nmrbox:
        try:
            cif_file_path = '/reboxitory/2021/07/PDB/data/structures/all/mmCIF/{}.cif.gz'.format(pdb_id.lower())
            pdb_data = get_coordinates(cif_file_path,auth_tag,nmrbox=nmrbox)
        except FileNotFoundError:
            pdb_data = None
    else:
        if not os.path.isdir('./data'):
            os.system('mkdir ./data')
        if not os.path.isdir('./data/PDB'):
            os.system('mkdir ./data/PDB')
        cif_file = './data/PDB/{}.cif'.format(pdb_id)
        if not os.path.isfile(cif_file):
            cmd = 'wget https://files.rcsb.org/download/{}.cif -O ./data/PDB/{}.cif'.format(pdb_id, pdb_id)
            os.system(cmd)
        pdb_data = get_coordinates('./data/PDB/{}.cif'.format(pdb_id),auth_tag)
    return pdb_data

def get_bmrb_data(bmrb_id: Union[str,int],auth_tag: Optional[bool] = False, nmrbox: Optional[bool] =False):
    """
    Extract the chemical shift data for a given BMRB ID
    :param bmrb_id: BMRB ID
    :param auth_tag: Use sequence information from Auth_seq_ID, Optional
    :param nmrbox: Instead of BMRB-APR, use the NMR-STAR file from NMRBOx reboxitory. Optional
    :return: tuple of amide chemical  shift as dictionary, aromatic chemical shift as dictionary, entity size
        and entity and entity assembly size
    """
    if nmrbox:
        try:
            str_file_path = '/reboxitory/2021/07/BMRB/macromolecules/bmr{}/bmr{}_3.str'.format(
                bmrb_id,bmrb_id
            )
            bmrb_data = get_chemical_shifts(str_file_path, auth_tag)
        except FileNotFoundError:
            bmrb_data = None
    else:
        if not os.path.isdir('./data'):
            os.system('mkdir ./data')
        if not os.path.isdir('./data/BMRB'):
            os.system('mkdir ./data/BMRB')
        str_file = './data/BMRB/{}.str'.format(bmrb_id)
        xx=0
        if not os.path.isfile(str_file):
            cmd = 'wget http://rest.bmrb.io/bmrb/{}/nmr-star3 -O ./data/BMRB/{}.str'.format(bmrb_id, bmrb_id)
            xx=os.system(cmd)
        if xx==0:
            try:
                bmrb_data = get_chemical_shifts('./data/BMRB/{}.str'.format(bmrb_id),auth_tag)
            except AttributeError:
                bmrb_data = None
        else:
            bmrb_data = None
    return bmrb_data

def get_z_score_full(res: str, amide_cs: float) -> float:
    """
    Calculates Z-score for a given amide chemical shift
    :param res: residue name
    :param amide_cs: chemical shift
    :return: Z-score
    """
    m = {'ALA':8.193,'ARG':8.242,'ASN':8.331,'ASP':8.300,'CYS':8.379,'GLN':8.216,'GLU':8.330,'GLY':8.327,
         'HIS':8.258,'ILE':8.263,'LEU':8.219,'LYS':8.175,'MET':8.258,'PHE':8.337,'SER':8.277,'THR':8.235,
         'TRP':8.270,'TYR':8.296,'VAL':8.273}
    sd = {'ALA':0.642,'ARG':1.064,'ASN':0.983,'ASP':0.592,'CYS':0.697,'GLN':0.657,'GLU':0.750,'GLY':0.770,
         'HIS':0.734,'ILE':0.694,'LEU':0.652,'LYS':0.670,'MET':1.277,'PHE':0.732,'SER':0.602,'THR':0.641,
         'TRP':0.782,'TYR':0.741,'VAL':0.795}
    try:
        z_score = (amide_cs - m[res]) / sd[res]
    except KeyError:
        z_score = 0.00
    return round(z_score,3)

def get_z_score(res: str, amide_cs: float) -> float:
    """
    Calculates Z-score for a given amide chemical shift
    :param res: residue name
    :param amide_cs: chemical shift
    :return: Z-score
    """
    m = {'ALA':8.194,'ARG':8.234,'ASN':8.324,'ASP':8.300,'CYS':8.386,'GLN':8.219,'GLU':8.330,'GLY':8.330,
         'HIS':8.247,'ILE':8.262,'LEU':8.215,'LYS':8.177,'MET':8.251,'PHE':8.335,'SER':8.277,'THR':8.232,
         'TRP':8.264,'TYR':8.289,'VAL':8.270}
    sd = {'ALA':0.577,'ARG':0.601,'ASN':0.610,'ASP':0.558,'CYS':0.670,'GLN':0.569,'GLU':0.576,'GLY':0.619,
         'HIS':0.666,'ILE':0.674,'LEU':0.627,'LYS':0.589,'MET':0.575,'PHE':0.710,'SER':0.568,'THR':0.610,
         'TRP':0.761,'TYR':0.721,'VAL':0.659}
    try:
        z_score = (amide_cs - m[res]) / sd[res]
    except KeyError:
        z_score = 0.00
    return round(z_score,3)

if __name__ == "__main__":
    print (get_bmrb_data('11086'))
    print (get_pdb_data('1WYO'))