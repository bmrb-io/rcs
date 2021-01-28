from protein_builder import *
from proteins import *
import requests
import os
from multiprocessing import Pipe, cpu_count
from os import _exit as child_exit
import time
import traceback


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



def make_entries_list(entries_dict):
    """Break dict of pdb_ids and associated bmrb_ids into list of pairs."""
    entries_list = []
    for pdb_id in entries_dict:
        for bmrb_id in entries_dict[pdb_id]:
            entries_list.append([pdb_id, bmrb_id])
    return entries_list

def add_to_proteins_dict(protein, proteins_dict):
    """Add protein to proteins_dict."""
    pdb_id = protein.pdb_id
    bmrb_id = protein.bmrb_id
    if pdb_id not in proteins_dict:
        proteins_dict[pdb_id] = {}
    proteins_dict[pdb_id][bmrb_id] = protein
    return proteins_dict

def add_to_exceptions_map(exception, pdb_id, bmrb_id, exceptions_map):
    """Add exception to exceptions_map."""
    if pdb_id not in exceptions_map:
        exceptions_map[pdb_id] = {}
    exceptions_map[pdb_id][bmrb_id] = exception
    return exceptions_map

def get_proteins_dict_multi(entries_dict, build_anyway=False):
    """Get proteins_dict using multiprocessing to enhance performance."""
    os.nice(19) # So that we don't take up too many resources
    entries_list = make_entries_list(entries_dict)
    proteins_dict = {}
    exceptions_map = {}

    processes = []
    num_threads = cpu_count()
    for thread in range(0, num_threads):
        # Set up the pipes
        parent_conn, child_conn = Pipe()
        # Start the process
        processes.append([parent_conn, child_conn])
        new_pid = os.fork()
        # Okay, we are the child
        if new_pid == 0:
            child_conn.send("ready")
            while True:
                parent_message = child_conn.recv()
                if parent_message == 'die':
                    child_conn.close()
                    parent_conn.close()
                    child_exit(0)
                # Do work based on parent_message
                pdb_id = parent_message[0]
                bmrb_id = parent_message[1]
                try:
                    protein = get_protein(pdb_id, bmrb_id, build_anyway)
                    # Tell our parent we are ready for the next job
                    child_conn.send([protein, pdb_id, bmrb_id])
                except KeyError:
                    exception = "PDB ID not found in RCSB"
                    child_conn.send([exception, pdb_id, bmrb_id])
                except AttributeError as err:
                    exception = "BMRB entry deprecated."
                    child_conn.send([exception, pdb_id, bmrb_id])
                except Exception as err:
                    print(pdb_id, bmrb_id, err)
                    child_conn.send([err, pdb_id, bmrb_id])
        # We are the parent, don't need the child connection
        else:
            child_conn.close()

    # Check if entries have completed by listening on the sockets
    while len(entries_list) > 0:
        time.sleep(0.001)
        # Poll for processes ready to listen
        for process in processes:
            if process[0].poll(): # if it has something to say, I think.
                data = process[0].recv() # results from above
                if data: #if data is not empty
                    if data != "ready":
                        protein = data[0]
                        pdb_id = data[1]
                        bmrb_id = data[2]
                        if isinstance(protein, Protein):
                            proteins_dict = add_to_proteins_dict(
                                protein, proteins_dict
                            )
                        else:
                            exceptions_map = add_to_exceptions_map(
                                protein, pdb_id, bmrb_id, exceptions_map
                            )
                process[0].send(entries_list.pop()) #sends pdb_id and bmrb_id I think
                break # force to reevaluate len(entries_list)

    # Reap the children
    for thread in range(0, num_threads):
        # Get the last ready message from the child
        data = processes[thread][0].recv()
        # Tell the child to shut down
        processes[thread][0].send("die")
        res = os.wait() 
        if data:
            protein = data[0]
            pdb_id = data[1]
            bmrb_id = data[2]
            if isinstance(protein, Protein):
                proteins_dict = add_to_proteins_dict(protein, proteins_dict)
            else:
                exceptions_map = add_to_exceptions_map(
                    protein, pdb_id, bmrb_id, exceptions_map
                )

    return proteins_dict, exceptions_map

entries_dict = get_all_entries()
proteins_dict, exceptions_map = get_proteins_dict_multi(entries_dict)
reasons_dict = {}
for pdb_id in exceptions_map:
    for bmrb_id in exceptions_map[pdb_id]:
        reason = exceptions_map[pdb_id][bmrb_id]
        if reason not in reasons_dict:
            reasons_dict[reason] = 0
        reasons_dict[reason] += 1

for reason in reasons_dict:
    print(reason, reasons_dict[reason])