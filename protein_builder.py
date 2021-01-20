from k_file_maker import *
from k_file_reader import *
from noes_builder import *
from ids_map import *

def build_protein(pdb_id, bmrb_id):
    """
    Create a Protein object with amides, aromatic protons, chemical shifts,
    and NOE restraints if possible. Otherwise return an appropriate exception.

    Keyword arguments:
    pdb_id -- str of the PDB ID
    bmrb_id -- str of the BMRB ID
    Returns
    protein -- either a Protein object or an exception raised by 
        make_protein_from_file() or add_restraints()
    """
    ring_current_object = RingCurrentEffect(pdb_id, bmrb_id)
    k_file_path = ring_current_object.calculate_ring_current_effects(
        pdb_id, bmrb_id
    )
    if "NO_" in k_file_path or 'NONE' in k_file_path:
        print(pdb_id, bmrb_id, " DIDN'T WORK")
    else:
        protein = make_protein_from_file(k_file_path)
        if isinstance(protein, Protein):
            protein = add_restraints(protein)
        else:
            print(pdb_id, bmrb_id, protein)

def build_all_proteins(map_filename):
    id_map = IDs_map(map_filename)
    for bmrb_id, pdb_id in id_map.ids_list[:100]:
        build_protein(pdb_id, bmrb_id)