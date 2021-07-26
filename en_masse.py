from protein_builder import *
from proteins import Protein
import requests
import pandas as pd
import os
from multiprocessing import Pipe, cpu_count
from os import _exit as child_exit
import time
import traceback
from typing import Dict, Tuple, List

def make_entries_list() -> List[Tuple[str, str]]:
    """Break dict of pdb_ids and associated bmrb_ids into list of pairs."""
    filename = os.path.join(
        '/reboxitory', '2021', '06', 'BMRB', 'relational_tables', 'nmr-star3.1', 
        'web.pdb_link.csv'
    )
    entries_df = pd.read_csv(filename)

    pdb_ids = list(entries_df['pdb_id'])
    bmrb_ids = list(entries_df['bmrb_id'])

    entries_list = list(zip(pdb_ids, bmrb_ids))

    return entries_list

def add_to_proteins_dict(
    protein: Protein, proteins_dict: Dict[str, Dict[str, Protein]]
) -> Dict[str, Dict[str, Protein]]:
    """Add protein to proteins_dict."""
    pdb_id = protein.pdb_id
    bmrb_id = protein.bmrb_id
    if pdb_id not in proteins_dict:
        proteins_dict[pdb_id] = {}
    proteins_dict[pdb_id][bmrb_id] = protein
    return proteins_dict

def add_to_exceptions_map(
    exception: str, pdb_id: str, bmrb_id: str, 
    exceptions_map: Dict[str, Dict[str, str]]
) -> Dict[str, Dict[str, str]]:
    """Add exception to exceptions_map."""
    if pdb_id not in exceptions_map:
        exceptions_map[pdb_id] = {}
    exceptions_map[pdb_id][bmrb_id] = exception
    return exceptions_map

def get_proteins_dict_multi(
    build_anyway: bool = False
) -> Tuple[Dict[str, Dict[str, Protein]], Dict[str, Dict[str, str]]]:
    """Get proteins_dict using multiprocessing to enhance performance."""
    os.nice(19) # So that we don't take up too many resources
    entries_list = make_entries_list()
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
                if pdb_id == 'e':
                    print("WHAT")
                try:
                    protein = get_protein(pdb_id, bmrb_id, build_anyway)
                    # Tell our parent we are ready for the next job
                    child_conn.send([protein, pdb_id, bmrb_id])

                except Exception as err:
                    err = traceback.format_exc()
                    err=str(err)
                    child_conn.send([err, pdb_id, bmrb_id])
        # We are the parent, don't need the child connection
        else:
            child_conn.close()

    # Check if entries have completed by listening on the sockets
    while len(entries_list) > 0:
        time.sleep(0.001)
        # Poll for processes ready to listen
        for i, process in enumerate(processes):
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
                            err = categorize_err(str(protein))
                            exceptions_map = add_to_exceptions_map(
                                err, pdb_id, bmrb_id, exceptions_map
                            )
                ids = entries_list.pop()
                process[0].send(ids) #sends pdb_id and bmrb_id I think
                break # force to reevaluate len(entries_list)

    # Reap the children
    for thread in range(0, num_threads):
        # Get the last ready message from the child
        data = processes[thread][0].recv()
        # Tell the child to shut down
        processes[thread][0].send("die")
        res = os.wait() 
        if data:
            if data != 'ready':
                protein = data[0]
                pdb_id = data[1]
                bmrb_id = data[2]
                if isinstance(protein, Protein):
                    proteins_dict = add_to_proteins_dict(protein, proteins_dict)
                else:
                    err = categorize_err(str(protein))
                    exceptions_map = add_to_exceptions_map(
                        err, pdb_id, bmrb_id, exceptions_map
                    )

    return proteins_dict, exceptions_map

def categorize_err(err):

    if 'No such file or directory' in err:
        if 'PDB' in err:
            if 'nmr_restraints_v2' in err:
                return 'Restraint file not in reboxitory'
            elif 'mmCIF' in err:
                return 'mmCIF file not in reboxitory'
            else:
                return 'Unknown file error'
        elif 'BMRB' in err:
            return 'STR file not in reboxitory'
        else:
            return 'Unknown file error'
    elif 'Permission denied' in err:
        return 'Permission denied error in reboxitory'
    else:
        return err
    
        

