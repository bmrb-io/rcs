import json

def find_ids_diff(filename_1, filename_2):

    filename_1 = f'{filename_1}_ids.json'
    filename_2 = f'{filename_2}_ids.json'

    with open(filename_1, 'r') as f1:
        ids_1 = json.load(f1)
    with open(filename_2, 'r') as f2:
        ids_2 = json.load(f2)
    
    diff = []
    for ids in ids_1:
        if ids not in ids_2:
            diff.append(ids)
    print(diff)

def find_pairs_diff(filename_1, filename_2):

    filename_1 = f'{filename_1}_pairs.json'
    filename_2 = f'{filename_2}_pairs.json'

    with open(filename_1, 'r') as f1:
        pairs_1 = json.load(f1)
    with open(filename_2, 'r') as f2:
        pairs_2 = json.load(f2)
    
    diff = []
    for pair in pairs_1:
        if pair not in pairs_2:
            diff.append(pair)
    print(diff)

    

find_pairs_diff('test_from_scratch', 'test_from_load')