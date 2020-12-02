from shift_values import *
from positions import *
from ids_map import *
import matplotlib.pyplot as plt

'''
shifts_list = get_shifts_list('4020')
atoms_dict = get_positions('1BRV')

pdb_aminos_list = list(atoms_dict.keys())
pdb_aminos_list = [(int(amino[0]), amino[1]) for amino in pdb_aminos_list]
pdb_aminos_list.sort()

bmrb_aminos_list = list(set([(int(shift[0]), shift[1]) for shift in shifts_list]))
bmrb_aminos_list.sort()

match_shifts(shifts_list, atoms_dict)
'''

ids_list = IDs_map("ids_list.json").ids_list
correctness_list = []
i = 0
list_len = len(ids_list)
for bmrb_id, pdb_id in ids_list[:200]:
    print(i)
    i += 1
    try:
        shifts_list = get_shifts_list(bmrb_id)
        if shifts_list is not None:
            atoms_dict = get_positions(pdb_id)
            nt, nf, framed_bool = match_shifts(shifts_list, atoms_dict)
            correctness_list.append(nt / (nt + nf))
    except OSError:
        print("Whoops! That entry isn't available in BMRB right now.", bmrb_id)

fig, ax = plt.subplots()
ax.hist(correctness_list, density=True, bins=30, rwidth=0.9)
plt.show()
