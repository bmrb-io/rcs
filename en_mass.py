from protein_builder import *
from proteins import *
import requests

def get_proteins_dict(entries_dict, build_anyway=False):
    """
    Create a Protein object for all viable entries, and return an exception 
    otherwise. Add successful proteins to proteins_dict and exceptions to
    exceptions_map.

    Keyword arguments:
    map_filename -- path to file containing map of bmrb_ids to pdb_ids
    Returns:
    proteins_dict -- dict of all successfully created proteins
    exceptions_map -- dict of reasons for failure to create proteins
    """
    proteins_dict = {}
    exceptions_map = {}
    for pdb_id in entries_dict:
        for bmrb_id in entries_dict[pdb_id]:
            print(pdb_id)
            protein = get_protein(pdb_id, bmrb_id, build_anyway)
            if isinstance(protein, Protein):
                if pdb_id not in proteins_dict:
                    proteins_dict[pdb_id] = {}
                proteins_dict[pdb_id][bmrb_id] = protein
            else:
                if pdb_id not in exceptions_map:
                    exceptions_map[pdb_id] = {}
                exceptions_map[pdb_id][bmrb_id] = protein
    return proteins_dict, exceptions_map

def get_all_entries():
    """Return a dict of all corresponding PDB and BMRB IDs."""
    url = "http://api.bmrb.io/v2/mappings/bmrb/pdb"
    r = requests.get(url).json()
    entries_dict = {}
    for ids_dict in r:
        bmrb_id = ids_dict['bmrb_id']
        pdb_ids = ids_dict['pdb_ids']
        for pdb_id in pdb_ids:
            if pdb_id not in entries_dict:
                entries_dict[pdb_id] = []
            entries_dict[pdb_id].append(bmrb_id)
    return entries_dict


import os
from multiprocessing import Pipe, cpu_count
from os import _exit as child_exit
import time

os.nice(19)

def get_proteins_dict_multi(entries_dict):
    processes = []
    num_threads = cpu_count()
    for thread in range(0, num_threads):
        # Set up the pipes
        parent_conn, child_conn = Pipe()
        # Start the process
        processes.append([parent_conn, child_conn])
        new_pid = os.fork()
        # Okay, we are the child
        if new_pid == 0;
            child_conn.send("ready")
            while True:
                parent_message = child_conn.recv()
                if parent_message == 'die':
                    child_conn.close()
                    parent_conn.close()
                    child_exit(0)
                # Do work based on parent_message
                result = get_protein(parent_message[0], parent_message[1])
                # Tell our parent we are ready for the next job
                child_conn.send(result)
        # We are the parent, don't need the child connection
        else:
            child_conn.close()

    # Check if entries have completed by listening on the sockets
    while len(to_process['combined']) > 0: ##WHERE DID THIS COME FROM
        time.sleep(0.001)
        # Poll for processes ready to listen
        for process in processes:
            if process[0].poll() ##WHAT IS THIS DOING
                data = process[0].recv()
                if data:
                    if data != "ready":  ##WHEN WOULD THIS HAPPEN
                        add_to_loaded(data) ##WHERE DID THIS COME FROM
                process[0].send(to_process['combined'].pop()) ##Come back when less sleepy
                break #in case this was the last child, I think

    # Reap the children
    for thread in range(0, num_threads):
        # Get the last ready message from the child
        data = processes[thread][0].recv()
        # Tell the child to shut down
        processes[thread][0].send("die")
        res = os.wait() ##to let everything catch up?
        if data:
            add_to_loaded(data)




