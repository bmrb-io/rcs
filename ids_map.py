import requests
import json


class IDs_map:

    def __init__(self, file_name, make_file = False):
        
        self.ids_list = []
        if make_file:
            self.make_ids_map_json(file_name)
        ids_map_json = self.get_ids_map_json(file_name)
        self.make_ids_list(ids_map_json)



    def make_ids_map_json(self, file_name):
        
        url = "http://api.bmrb.io/v2/mappings/bmrb/pdb"
        r = requests.get(url).json()
        with open(file_name, 'w') as f:
            json.dump(r, f)
        

    def get_ids_map_json(self, file_name):

        with open(file_name) as f:
            ids_map_json = json.load(f)
        return ids_map_json
    
    def make_ids_list(self, ids_map_json):

        for item in ids_map_json:
            bmrb_id = item['bmrb_id']
            pdb_id = item['pdb_ids'][0] #for now, we only consider one
            self.ids_list.append((bmrb_id, pdb_id))


        