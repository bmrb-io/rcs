from k_file_maker import *
from k_file_reader import *
from noes_builder import *

def build_protein(pdb_id, bmrb_id):
    ring_current_object = RingCurrentEffect(pdb_id, bmrb_id)
    k_file_path = ring_current_object.calculate_ring_current_effects(
        pdb_id, bmrb_id
    )
    protein = make_protein_from_file(k_file_path)
    if isinstance(protein, Protein):
        protein = add_restraints(protein)
    return protein

#build_protein('2L4N', '17245')