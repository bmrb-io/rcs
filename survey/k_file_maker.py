from pdbecif.mmcif_io import CifFileReader
from math import sqrt, acos
import numpy
import operator
import pynmrstar
from numpy import mean
import os,sys
import csv
import os.path
from operator import itemgetter
import requests

class RingCurrentEffect(object):
    atoms = {
        'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', 'HD1', 'HD2', 'HE1', 'HE2', 'xx', 'yy']  # if needed un comment
    }


    def __init__(self,pdbid,bmrbid):
        #self.cal_prop(bmrbid)
        # self.calculate_ring_current_effects(pdbid,bmrbid)
        #self.generate_job_files(1000)
        pass

    def cal_mean_distance(self,pdb,atm1,atm2):
        d=[]
        for m in pdb.keys():
            try:
                d.append(self.get_distance(pdb[m][atm1],pdb[m][atm2]))
            except KeyError:
                d.append(0.0)
        return len(d),numpy.mean(d),numpy.std(d)

    def get_seq(self,str_file):
        str_data = pynmrstar.Entry.from_file(str_file)
        sq_dat = str_data.get_tag('_Entity_comp_index.Comp_ID')
        sqq_dat = str_data.get_tag('_Entity.Polymer_type')
        n=len(sq_dat)

        return n,sqq_dat,[sq_dat.count('HIS'),sq_dat.count('TYR'),sq_dat.count('PHE'),sq_dat.count('TRP')]




    def calculate_ring_current_effects(self,pdbid,bmrbid):
        '''
        Extract the information about amide chemical shift and the nearby aromatic ring geometry and chemical shifts and write out csv file
        :param pdbid: matching pdb id example 2L4N
        :param bmrbid: matching bmrb id example 17245
        '''

        cif_file = os.path.join(
            '/reboxitory', '2021', '06', 'PDB', 'data', 'structures', 'all', 
            'mmCIF', f'{pdbid.lower()}.cif.gz'
        )
        str_file = os.path.join(
            '/reboxitory', '2021', '06', 'BMRB', 'macromolecules',
            f'bmr{bmrbid}', f'bmr{bmrbid}_3.str'
        )

        if not os.path.isfile(cif_file):
            raise ValueError('mmCIF file not in reboxitory')
        if not os.path.isfile(str_file):
            raise ValueError('STR file not in reboxitory')
        pdb_auth = self.get_coordinates(cif_file,pdbid,use_auth_tag=True) # creates the dictionary using original seq no
        pdb_orig = self.get_coordinates(cif_file,pdbid,use_auth_tag=False) # creates the dictionary using author seq no
        pdb_auth_keys = [i for i in pdb_auth[0][1].keys() if i[3]=='H']
        pdb_orig_keys = [i for i in pdb_orig[0][1].keys() if i[3]=='H']
        cs_orig = self.get_chemical_shifts(str_file,auth_tag=False) # creates the seq using original seq no
        cs_auth = self.get_chemical_shifts(str_file,auth_tag=True) # creates the seq using original seq no
        cs_orig_keys = [i for i in cs_orig[0].keys() if i[3]=='H']
        cs_auth_keys = [i for i in cs_auth[0].keys() if i[3] == 'H']
        orig_auth = [i for i in cs_orig_keys if i in pdb_auth_keys]
        orig_orig = [i for i in cs_orig_keys if i in pdb_orig_keys]
        auth_auth = [i for i in cs_auth_keys if i in pdb_auth_keys]
        auth_orig = [i for i in cs_auth_keys if i in pdb_orig_keys]
        pdb = pdb_orig
        if len(cs_orig_keys)>0:
            seq_len=float(len(cs_orig_keys))
            orig_orig_match=float(len(orig_orig))/seq_len
            orig_auth_match=float(len(orig_auth))/seq_len
            auth_auth_match=float(len(auth_auth))/seq_len
            auth_orig_match=float(len(auth_orig))/seq_len

            if orig_orig_match > 0.7:
                pdb = pdb_orig
                cs = cs_orig
                tag_match='ORIG_ORIG'
            elif orig_auth_match > 0.7:
                pdb = pdb_auth
                cs = cs_orig
                tag_match = 'ORIG_AUTH'
            elif auth_auth_match > 0.7:
                pdb = pdb_auth
                cs = cs_auth
                tag_match = 'AUTH_AUTH'
            elif auth_orig_match > 0.7:
                pdb = pdb_orig
                cs = cs_auth
                tag_match = 'AUTH_ORIG'
            else:
                pdb = None
                cs = cs_orig
                tag_match = None
        else:
            pdb = None
            cs = [None,None,None,None]
            tag_match = None

        if tag_match is not None:
            ar=self.find_aromatic_residues(pdb[0])
        else:
            ar = None

        fout = self.find_amide_ring_distance(pdb, ar, cs , pdbid, bmrbid, tag_match)
        return fout




    @staticmethod
    def get_distance(c1, c2):
        """
        Calculates the distance between two coordinate points
        :param c1: array of x,y,z
        :param c2: array of x,y,z
        :return: distance between two ponts
        """
        return numpy.linalg.norm(c1 - c2)

    @staticmethod
    def get_chemical_shifts(str_file,auth_tag=False):
        atoms_cs = {
            'PHE': ['CG','CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
            'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
            'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1'],
            'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', 'HD1', 'HD2', 'HE1', 'HE2', 'xx', 'yy']  # if needed un comment
        }
        str_data = pynmrstar.Entry.from_file(str_file)
        csdata = str_data.get_loops_by_category('Atom_chem_shift')
        assembly_data = str_data.get_loops_by_category('Entity_assembly')
        entity = str_data.get_tag('_Entity_assembly.Entity_ID')
        entity_size=len(set(entity))
        assembly_size=(len(entity))
        csh = {}
        csh2 = {}
        for cs in csdata:
            tag_list = cs.get_tag_names()
            id2 = tag_list.index('_Atom_chem_shift.Auth_asym_ID')
            if auth_tag:
                id1 = tag_list.index('_Atom_chem_shift.Auth_seq_ID')
                id3 = tag_list.index('_Atom_chem_shift.Auth_comp_ID')
            else:
                id1 = tag_list.index('_Atom_chem_shift.Comp_index_ID')
                id3 = tag_list.index('_Atom_chem_shift.Comp_ID')
            id3 = tag_list.index('_Atom_chem_shift.Comp_ID')
            id4 = tag_list.index('_Atom_chem_shift.Atom_ID')
            id5 = tag_list.index('_Atom_chem_shift.Val')
            id6 = tag_list.index('_Atom_chem_shift.Ambiguity_code')
            for d in cs.data:
                if d[id4] == 'H':
                    if d[id2] in ['.','1']:
                        d[id2] = 'A' # temp fix
                    csh[(d[id1],d[id2],d[id3],d[id4])]=d[id5]
                if d[id3] in atoms_cs.keys():
                    if d[id2] in ['.','1']:
                        d[id2] = 'A' #temp fix
                    k = (d[id1],d[id2],d[id3])
                    if d[id4] in atoms_cs[d[id3]]:
                        if k not in csh2.keys():
                            csh2[k] = {}
                        csh2[k][d[id4]]=(d[id5],d[id6])
        return csh,csh2,entity_size,assembly_size

    @staticmethod
    def solid_angle(a_deg,r):
        s=1.4
        #A=((3.0*numpy.sqrt(3))/2.0)*s*s
        a=(numpy.pi/180)*a_deg
        r1=r*1e10
        #sa2=2*numpy.pi*(1.0-1.0/(numpy.sqrt(1+(A*numpy.cos(a)/(numpy.pi*r1**r1)))))
        sa = 2 * numpy.pi * (1.0 - 1.0 / (numpy.sqrt(1 + ( numpy.cos(a) / (numpy.pi * r ** r)))))
        #print (a_deg)
        sa_deg=(180.0/numpy.pi)*sa
        #sa_deg2 = (180.0 / numpy.pi) * sa2
        #print (a_deg,sa_deg,sa_deg2)
        return sa_deg

    @staticmethod
    def get_coordinates(cif_file, pdbid, use_auth_tag=True):
        """
        Extract coordinate information from cif file as a dictionary
        {model_id : {(seq_id,chain_id,res_id,atom_id) : array[x,y,x],...},...}
        :param cif_file: Input coordinate file
        :return: dictionary
        """
        cfr = CifFileReader()
        cif_obj = cfr.read(cif_file, output='cif_wrapper')
        cif_data = cif_obj[pdbid]
        atom_sites = cif_data._atom_site
        model_nums = set(atom_sites['pdbx_PDB_model_num'])
        pdb_models = {}
        atom_ids = {}
        for i in model_nums:
            pdb_models[int(i)] = {}
            atom_ids[int(i)] = {}
        for atom_site in list(atom_sites):
            model_id = int(atom_site['pdbx_PDB_model_num'])
            aut_seq_id = atom_site['auth_seq_id']
            aut_asym_id = atom_site['auth_asym_id']
            aut_comp_id = atom_site['auth_comp_id']
            aut_atom_id = atom_site['auth_atom_id']
            entity_id = atom_site['label_entity_id']
            asym_id = atom_site['label_asym_id']
            seq_id = atom_site['label_seq_id']
            comp_id = atom_site['label_comp_id']
            atom_id = atom_site['label_atom_id']
            alt_id = atom_site['label_alt_id']
            icode_id = atom_site['pdbx_PDB_ins_code']
            posn_x = atom_site['Cartn_x']
            posn_y = atom_site['Cartn_y']
            posn_z = atom_site['Cartn_z']
            posn = numpy.array([float(posn_x), float(posn_y), float(posn_z)])
            if use_auth_tag:
                key = (aut_seq_id, aut_asym_id, aut_comp_id, aut_atom_id)
            else:
                key = (seq_id, asym_id, comp_id, atom_id)
            val_atom = (
                entity_id, asym_id, comp_id, seq_id, aut_seq_id, alt_id,
                icode_id, aut_asym_id
            )
            val_pdb = posn
            atom_ids[model_id][key] = val_atom
            pdb_models[model_id][key] = val_pdb
        return pdb_models, atom_ids

    @staticmethod
    def get_sigma_value(res, x):
        m = {
            'ALA': 8.194, 'ARG': 8.234, 'ASN': 8.324, 'ASP': 8.300, 'CYS': 8.386,
            'GLN': 8.219, 'GLU': 8.330, 'GLY': 8.330, 'HIS': 8.247, 'ILE': 8.262,
            'LEU': 8.215, 'LYS': 8.177, 'MET': 8.251, 'PHE': 8.335, 'SER': 8.277,
            'THR': 8.232, 'TRP': 8.264, 'TYR': 8.289, 'VAL': 8.270,
        }
        sd = {
            'ALA': 0.577, 'ARG': 0.601, 'ASN': 0.610, 'ASP': 0.558, 'CYS': 0.670,
            'GLN': 0.569, 'GLU': 0.576, 'GLY': 0.619, 'HIS': 0.666, 'ILE': 0.674,
            'LEU': 0.627, 'LYS': 0.589, 'MET': 0.575, 'PHE': 0.710, 'SER': 0.568,
            'THR': 0.610, 'TRP': 0.761, 'TYR': 0.721, 'VAL': 0.659
        }
        try:
            sp = (x-m[res])/sd[res]
        except KeyError:
            sp = 0.00
        return round(sp,3)

    @staticmethod
    def get_centroid(p):
        #print (len(p),p)
        x= [i[0] for i in p]
        y= [i[1] for i in p]
        z= [i[2] for i in p]
        c = [sum(x)/len(p),sum(y)/len(p),sum(z)/len(p)]
        return numpy.array(c)

    def find_angle(self,p,c,cn,d):
        pc = self.get_centroid(p)
        v1=p[1]-pc
        v2=p[2]-pc
        nv=numpy.cross(v1,v2)
        nv2 = numpy.cross(v2,v1)
        cv=c-pc
        cv2=c-cn
        nnv=nv/numpy.linalg.norm(nv)
        ncv=cv/numpy.linalg.norm(cv)
        ncv2 = cv2 / numpy.linalg.norm(cv2)
        dp = abs(numpy.dot(nnv,ncv))
        dp2 = abs(numpy.dot(nnv, ncv2))
        ang = numpy.arccos(dp)
        ang2 = numpy.arccos(dp2)
        ang_deg = (180/numpy.pi)*ang
        ang_deg2 = (180 / numpy.pi) * ang2
        #print(ang_deg)
        s_ang=self.solid_angle(ang_deg,d*1e-10)
        return ang_deg,s_ang



    def find_mean_distance(self,ph,pc):
        d = []
        for d1 in pc:
            d.append(numpy.linalg.norm(ph - d1))
        return mean(d), numpy.std(d)


    def find_aromatic_residues(self,pdb):
        rl = []
        for a in pdb[1].keys():
            if a[2] in ['PHE','TRP','TYR','HIS']:
                rl.append((a[0],a[1],a[2]))
        return list(set(rl))



    def find_amide_ring_distance(self,pdb2,aromatic,cs2,pdbid,bmrbid,tag_match):
        atoms_cs = {
            'PHE': ['CG','CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
            'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
            'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1'],
           'HIS' :['CG','ND1','CD2','CE1','NE2','HD1','HD2','HE1','HE2','xx','yy'] #if needed un comment
        }
        amide_chemical_shift,aromatic_chemical_shift,entity_size,assembly_size = cs2
        if not os.path.isdir('./output'):
            os.system('mkdir ./output')
        if tag_match is not None:
            fout = './output/{}_{}_{}.dat'.format(pdbid, bmrbid, tag_match)
        else:
            fout = './output/{}_{}_NONE.dat'.format(pdbid, bmrbid) #seq mismatch case
        with open(fout, 'w') as fo:
            if tag_match is None:
                fo.write('{},{} Severe residue index mismatch'.format(pdbid, bmrbid))
            else:
                pdb = pdb2[0]
                aromatic_atoms = {}
                for a in aromatic:
                    aromatic_atoms[a] = []
                    for a1 in self.atoms[a[2]]:
                        aromatic_atoms[a].append((a[0], a[1], a[2], a1))
                if len(aromatic) == 0:
                    fo.write('{},{} No aromatic residues in sequence'.format(pdbid, bmrbid)) #Not deposited or DNA/RNA?
                elif len(amide_chemical_shift)==0:
                    fo.write('{},{} No amide shifts reported'.format(pdbid, bmrbid)) #DNA/RNA and other non-polypepdie
                else:
                    for a in pdb[1].keys():
                        if a[3] == 'H':
                            ar_info=[]
                            for ar in aromatic_atoms.keys():
                                md=[]
                                sd=[]
                                cd = []
                                ang = []
                                ang2 = []
                                for m in pdb.keys():
                                    p = []
                                    patom=[]
                                    for atom in aromatic_atoms[ar]:
                                        try:
                                            patom.append(atom)
                                            p.append(pdb[m][atom])
                                        except KeyError:
                                            pass
                                    c = self.get_centroid(p)
                                    d = self.get_distance(pdb[m][a],c)
                                    mean_d,std_d = self.find_mean_distance(pdb[m][a],p)
                                    md.append(mean_d)
                                    sd.append(std_d)
                                    cd.append(d)
                                    an=(a[0],a[1],a[2],'N')
                                    #print (a,ar,m)
                                    angles = self.find_angle(p, pdb[m][a],pdb[m][an],d)
                                    ang.append(angles[0])
                                    ang2.append(angles[1])
                                ar_info.append([round(numpy.mean(cd),3),
                                            round(numpy.std(cd),3),
                                            round(numpy.mean(ang),3),
                                            round(numpy.std(ang),3),
                                            round(numpy.mean(ang2), 3),
                                            round(numpy.std(ang2), 3),ar])
                            arr_info=sorted(ar_info, key=lambda x: x[0])
                            m = {'ALA': 8.193, 'ARG': 8.242, 'ASN': 8.331, 'ASP': 8.300, 'CYS': 8.379, 'GLN': 8.216,
                                'GLU': 8.330, 'GLY': 8.327,
                                'HIS': 8.258, 'ILE': 8.263, 'LEU': 8.219, 'LYS': 8.175, 'MET': 8.258, 'PHE': 8.337,
                                'SER': 8.277, 'THR': 8.235,
                                'TRP': 8.270, 'TYR': 8.296, 'VAL': 8.273}
                            if a in amide_chemical_shift.keys() and a[2] in m.keys():
                                outdat='{},{},{},{},{},{},{},{}'.format(pdbid,bmrbid,a[0],a[2],amide_chemical_shift[a],self.get_sigma_value(a[2],float(amide_chemical_shift[a])),entity_size,
                                                                                                            assembly_size)
                                for kk in range(100):

                                    try:
                                        odat='{},{},{},{},{},{},{},{}'.format(arr_info[kk][-1][0],arr_info[kk][-1][2],
                                                                            arr_info[kk][0],arr_info[kk][1],
                                                                            arr_info[kk][2],arr_info[kk][3],
                                                                            arr_info[kk][4],arr_info[kk][5],
                                                                            arr_info[kk][6])
                                        close_atom_cs2 = []
                                        close_atom2=arr_info[kk][-1]
                                        for atm in atoms_cs[close_atom2[2]]:
                                            d_info=self.cal_mean_distance(pdb, a, (arr_info[kk][-1][0],arr_info[kk][-1][1],arr_info[kk][-1][2],atm))
                                            try:
                                                close_atom_cs2.append(aromatic_chemical_shift[close_atom2][atm][0])
                                                close_atom_cs2.append(aromatic_chemical_shift[close_atom2][atm][1])
                                            except KeyError:
                                                close_atom_cs2.append(".")
                                                close_atom_cs2.append(".")
                                            close_atom_cs2.append('{}'.format(round(d_info[1],3)))
                                            close_atom_cs2.append('{}'.format(round(d_info[2], 3)))
                                        close_atom_cs_all2 = ",".join(close_atom_cs2)
                                        odat='{},{}'.format(odat,close_atom_cs_all2)
                                    except IndexError:
                                        odat='.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.'
                                    outdat='{},{}'.format(outdat,odat)
                                outdat='{}\n'.format(outdat)
                                fo.write(outdat)
        return fout

    @staticmethod
    def generate_job_files(n):
        """Return a dict of all corresponding PDB and BMRB IDs."""
        url = "http://api.bmrb.io/v2/mappings/bmrb/pdb?match_type=exact"
        r = requests.get(url).json()
        i=0
        j=0
        for ids_dict in r:
            if j==n:
                j=0
                f.close()
            if j==0 :
                i += 1
                f=open('job_{}.sh'.format(i),'w')
                f.write('#!/usr/bin/env bash\n')

            bmrb_id = ids_dict['bmrb_id']
            pdb_ids = ids_dict['pdb_ids']
            for k in pdb_ids:
                f.write('python3 RingCurrentEffects.py {} {}\n'.format(bmrb_id,k))
                j+=1


    @staticmethod
    def get_pdb(pdb_id):
        cmd = 'wget https://files.rcsb.org/download/{}.cif -O ./data/PDB/{}.cif'.format(pdb_id,pdb_id)
        os.system(cmd)

    @staticmethod
    def get_bmrb(bmrb_id):
        cmd = 'wget http://rest.bmrb.io/bmrb/{}/nmr-star3 -O ./data/BMRB/{}.str'.format(bmrb_id,bmrb_id)
        os.system(cmd)

    @staticmethod
    def check_output(datafile):
        id_pair=[]
        url = "http://api.bmrb.io/v2/mappings/bmrb/pdb?match_type=exact"
        r = requests.get(url).json()
        for ids_dict in r:
            bmrb_id = ids_dict['bmrb_id']
            pdb_ids = ids_dict['pdb_ids']
            for k in pdb_ids:
                id_pair.append('{}-{}'.format(bmrb_id,k))
        f = open('missing_job.sh', 'w')
        f.write('#!/usr/bin/env bash\n')
        #print (len(id_pair))
        id_pair2=[]
        with open(datafile) as csvfile:
            readcsv = csv.reader(csvfile, delimiter=',')
            for d in readcsv:
                if len(d[1]) > 6:
                    bid=d[1].split(" ")[0]
                else:
                    bid = d[1]
                id='{}-{}'.format(bid,d[0])

                id_pair2.append(id)
                #print (id,id_pair.index(id))
                    #print (id)
        u_ids=list(set(id_pair2))
        n_ids=[]
        for id in id_pair:
            if id not in u_ids:
                f.write('python3 RingCurrentEffects.py {} {}\n'.format(id.split("-")[0], id.split("-")[1]))
                #print (id)
                n_ids.append(id)
        f.close()
        #print (len(id_pair))
        #print (len(u_ids))
        #print (len(n_ids))


'''
if __name__ == "__main__":
    bmrbid = sys.argv[1]
    pdbid = sys.argv[2]
    # bmrbid = '11086'
    # pdbid = '5KVP'
    # bmrbid = '30139'
    # pdbid = '2JO7'
    # pdbid= '2KAH'
    # bmrbid='16023'
    p=RingCurrentEffect(pdbid,bmrbid)
    #p.check_output('data_08022021_3.csv')
'''
