import requests
import json
from aminos import Amino
from atoms import Atom

class PDB_map:

    def __init__(self, file_name, make_file = False):
        self.file_name = file_name
        self.ids_dict = {}
        self.pdb_map_json = None

        if make_file:
            self.make_pdb_map_json()
        self.get_pdb_map_json()
        self.make_ids_dict()



    def make_pdb_map_json(self):
        
        url = "http://api.bmrb.io/v2/mappings/bmrb/pdb"
        r = requests.get(url).json()
        with open(self.file_name, 'w') as f:
            json.dump(r, f)
        

    def get_pdb_map_json(self):

        with open(self.file_name) as f:
            self.pdb_map_json = json.load(f)
    
    def make_ids_dict(self):

        for item in self.pdb_map_json:
            bmrb_id = item['bmrb_id']
            pdb_id = item['pdb_ids'][0] #for now, we only consider one
            self.ids_dict[bmrb_id] = pdb_id


        