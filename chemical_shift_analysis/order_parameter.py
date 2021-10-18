import csv
import numpy
import json
import sys
from urllib.request import urlopen, Request

_API_URL = "http://api.bmrb.io/v2"


def get_order_param(eid):
    url = Request(_API_URL + "/entry/{}?loop=Order_param".format(eid))
    print (_API_URL + "/entry/{}?loop=Order_param".format(eid))
    url.add_header('Application', 'PyBMRB')
    r = urlopen(url)
    dump = json.loads(r.read())
    #print (dump)
    od={}
    if len(dump[eid]['Order_param']) == 0:
        print ('No order parameter for {}'.format(eid))
    else:
        for d in dump[eid]['Order_param']:
            seq_id = d['tags'].index('Comp_index_ID')
            comp_id = d['tags'].index('Comp_ID')
            atm_id = d['tags'].index('Atom_ID')
            od_pm=d['tags'].index('Order_param_val')
            for i in d['data']:
                k='{}-{}-{}'.format(i[seq_id],i[comp_id],i[atm_id])
                if k not in od.keys():
                    od[k]=[]
                od[k].append(i[od_pm])
    return od


def append_order_parameter(datafile):
    with open(datafile) as csvfile:
        mean_d=[]
        std_d=[]
        z=[]
        cs=[]
        mean_ang=[]
        std_ang=[]
        lab=[]
        color=[]
        sigma=[]
        ccc=[]
        sigma_cutoff=3.0
        pdb_total = []
        pdb_subset =[]
        solid_ang=[]
        sa=[]
        ret=[]
        readcsv = csv.reader(csvfile, delimiter=',')
        i=0
        od='999999'
        fo=open('data_woth_odpm_H.csv','w')
        for data in readcsv:
            amide_info=data[0:8]
            eid = amide_info[1]
            if eid != 'bmrb_id':
                k='{}-{}-{}'.format(amide_info[2],amide_info[3],'H')
                if od != eid:
                    ord_param = get_order_param(eid)
                    if len(ord_param)>0:
                        print (ord_param)
                    od = eid
                if len(ord_param)>0:
                    try:
                        odpm=ord_param[k][0]
                    except KeyError:
                        odpm='.'
                else:
                    odpm='.'
                #print (odpm)
                data.append(odpm)
                fo.write('{}\n'.format(",".join(data)))
        fo.close()

if __name__ == "__main__":
    csvfile = sys.argv[1]
    append_order_parameter(csvfile)

