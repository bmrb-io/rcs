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

class RingCurrentEffect(object):
    atoms = {
        'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', 'HD1', 'HD2', 'HE1', 'HE2', 'xx', 'yy']  # if needed un comment
    }

    def __init__(self,pdbid,bmrbid):
        pass

    def cal_mean_distance(self,pdb,atm1,atm2):
        d=[]
        for m in pdb.keys():
            try:
                d.append(self.get_distance(pdb[m][atm1],pdb[m][atm2]))
            except KeyError:
                d.append(0.0)
        return len(d),numpy.mean(d),numpy.std(d)

    def cal_prop(self,bmrbid):
        print ('Calculating aromatic ring and amide proton interaction in {} {}'.format(pdbid, bmrbid))
        if not os.path.isdir('./data'):
            os.system('mkdir ./data')
        if not os.path.isdir('./data/BMRB'):
            os.system('mkdir ./data/BMRB')
        cif_file = './data/PDB/{}.cif'.format(pdbid)
        str_file = './data/BMRB/{}.str'.format(bmrbid)
        if not os.path.isfile(cif_file):
            self.get_pdb(pdbid)
        if not os.path.isfile(str_file):
            self.get_bmrb(bmrbid)
        n,sq,x=self.get_seq(str_file)
        if not os.path.isdir('./output'):
            os.system('mkdir ./output')
        fout = './output/{}_{}.sq'.format(pdbid, bmrbid)
        fo = open(fout, 'w')
        fo.write("{},{},{},{},{},{},{},{},{}\n".format(bmrbid,pdbid,x[0],x[1],x[2],x[3],sum(x),n,",".join(sq)))
        fo.close()

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
        if not os.path.isdir('./data'):
            os.system('mkdir ./data')
        if not os.path.isdir('./data/PDB'):
            os.system('mkdir ./data/PDB')
        if not os.path.isdir('./data/BMRB'):
            os.system('mkdir ./data/BMRB')
        cif_file = './data/PDB/{}.cif'.format(pdbid)
        if not os.path.isfile(cif_file):
            self.get_pdb(pdbid)
        pdb_auth = self.get_coordinates(cif_file, pdbid, use_auth_tag=True) # creates the dictionary using original seq no
        pdb_orig = self.get_coordinates(cif_file, pdbid, use_auth_tag=False) # creates the dictionary using author seq no
        auth_keys = [i for i in pdb_auth[0][1].keys() if i[3]=='H']
        orig_keys = [i for i in pdb_orig[0][1].keys() if i[3]=='H']
        cs = self.get_chemical_shifts(bmrbid) # creates the seq using original seq no
        cs_keys = [i for i in cs[0].keys() if i[3]=='H']
        key_match_auth = [i for i in cs_keys if i in auth_keys]
        key_match_orig = [i for i in cs_keys if i in orig_keys]
        pdb = pdb_orig
        if len(auth_keys)==0:
            pdb = None
        else:
            auth_math = float(len(key_match_auth)) / float(len(auth_keys))
        if len(orig_keys)==0:
            pdb = None
        else:
            orig_match = float(len(key_match_orig)) / float(len(orig_keys))
        # check the seq numbering in pdb and bmrb matches at lease 70% and identify which tag is matching
        if pdb is not None:
            if orig_match > 0.7:
                tag_match = 'ORIG'
                pdb = pdb_orig
            elif auth_math > 0.7:
                tag_match = 'AUTH'
                pdb = pdb_auth
            else:
                tag_match = None
                pdb = None
            if tag_match is not None:
                ar=self.find_aromatic_residues(pdb[0])
            else:
                ar = None
        else:
            ar=None
            tag_match = None
        fout = self.find_amide_ring_distance(pdb, ar, cs, pdbid, bmrbid, tag_match)
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
    def get_chemical_shifts(bmrbid):
        atoms_cs = {
            'PHE': ['CG','CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
            'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
            'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1'],
            'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', 'HD1', 'HD2', 'HE1', 'HE2', 'xx', 'yy']  # if needed un comment
        }
        str_file = f'./data/BMRB/{bmrbid}.str'
        if not os.path.isfile(str_file):
            str_data = pynmrstar.Entry.from_database(bmrbid)
            str_data.write_to_file(str_file)
        else:
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
            id1 = tag_list.index('_Atom_chem_shift.Comp_index_ID')
            id2 = tag_list.index('_Atom_chem_shift.Auth_asym_ID')
            id3 = tag_list.index('_Atom_chem_shift.Comp_ID')
            id4 = tag_list.index('_Atom_chem_shift.Atom_ID')
            id5 = tag_list.index('_Atom_chem_shift.Val')
            id6 = tag_list.index('_Atom_chem_shift.Ambiguity_code')
            for d in cs.data:
                if d[id4] == 'H':
                    if d[id2] == '.':
                        d[id2] = 'A' # temp fix
                    csh[(d[id1],d[id2],d[id3],d[id4])]=d[id5]
                if d[id3] in atoms_cs.keys():
                    if d[id2] == '.':
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
        A=((3.0*numpy.sqrt(3))/2.0)*s*s
        a=(numpy.pi/180)*a_deg
        r1=r*1e10
        sa2=2*numpy.pi*(1.0-1.0/(numpy.sqrt(1+(A*numpy.cos(a)/(numpy.pi*r1**r1)))))
        sa = 2 * numpy.pi * (1.0 - 1.0 / (numpy.sqrt(1 + ( numpy.cos(a) / (numpy.pi * r ** r)))))
        #print (a_deg)
        sa_deg=(180.0/numpy.pi)*sa
        sa_deg2 = (180.0 / numpy.pi) * sa2
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
        max_models = int(atom_sites['pdbx_PDB_model_num'][-1])
        pdb_models = {}
        atom_ids = {}
        for i in range(max_models):
            pdb_models[i+1] = {}
            atom_ids[i+1] = {}
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
        m = {'ALA':8.193,'ARG':8.242,'ASN':8.331,'ASP':8.300,'CYS':8.379,'GLN':8.216,'GLU':8.330,'GLY':8.327,
             'HIS':8.258,'ILE':8.263,'LEU':8.219,'LYS':8.175,'MET':8.258,'PHE':8.337,'SER':8.277,'THR':8.235,
             'TRP':8.270,'TYR':8.296,'VAL':8.273}
        sd = {'ALA':0.642,'ARG':1.064,'ASN':0.983,'ASP':0.592,'CYS':0.697,'GLN':0.657,'GLU':0.750,'GLY':0.770,
             'HIS':0.734,'ILE':0.694,'LEU':0.652,'LYS':0.670,'MET':1.277,'PHE':0.732,'SER':0.602,'THR':0.641,
             'TRP':0.782,'TYR':0.741,'VAL':0.795}
        try:
            sp = (x-m[res])/sd[res]
        except KeyError:
            sp = 0.00
        return round(sp,3)

    @staticmethod
    def get_centroid(p):
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
        fo = open(fout, 'w')
        if tag_match is None:
            fo.write('{},{} No matching atoms found'.format(pdbid, bmrbid))
        else:
            pdb = pdb2[0]
            aromatic_atoms = {}
            for a in aromatic:
                aromatic_atoms[a] = []
                for a1 in self.atoms[a[2]]:
                    aromatic_atoms[a].append((a[0], a[1], a[2], a1))
            if len(aromatic) == 0:
                fo.write('{},{} No aromatic residues found'.format(pdbid, bmrbid)) #Not deposited or DNA/RNA?
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
                            for kk in range(5):
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
        fo.close()
        return fout
    @staticmethod
    def get_pdb(pdb_id):
        cmd = 'wget https://files.rcsb.org/download/{}.cif -O ./data/PDB/{}.cif'.format(pdb_id,pdb_id)
        os.system(cmd)

    @staticmethod
    def get_bmrb(bmrb_id):
        cmd = 'wget http://rest.bmrb.io/bmrb/{}/nmr-star3 -O ./data/BMRB/{}.str'.format(bmrb_id,bmrb_id)