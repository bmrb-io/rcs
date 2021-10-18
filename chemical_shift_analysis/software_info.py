import sys
import plotly.express as px
import pynmrstar
import csv
# import matplotlib.pyplot as plt
# from plotly.subplots import make_subplots
# import plotly.graph_objects as go
import numpy
def plot_software_info(datafile):
    def get_software_info(bmrbid):
        ent = pynmrstar.Entry.from_database(bmrbid)
        return ent.get_tag('_Software.Name')

    c1_list = ["rgb(0,206,227)", "rgb(20,206,227)", "rgb(40,206,227)", "rgb(60,206,227)", "rgb(80,206,227)",
               "rgb(100,206,227)"]
    c2_list = ["rgb(277,0,227)", "rgb(277,20,227)", "rgb(277,40,227)", "rgb(277,60,227)", "rgb(277,80,227)",
               "rgb(277,100,227)"]
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
        sigma_cutoff=2.0
        pdb_total = []
        pdb_subset =[]
        solid_ang=[]
        sa=[]
        ret=[]
        sa_std=[]
        d_cutoff=[]
        sw_info=[]
        bid=None
        readcsv = csv.reader(csvfile, delimiter=',')
        i=0
        for data in readcsv:
            amide_info=data[0:8]
            aro1_info=data[8:38]
            aro2_info = data[38:68]
            aro3_info = data[68:98]
            aro4_info = data[98:128]
            aro5_info = data[128:158]
            pdb_total.append(amide_info[0])
            if amide_info[1]!= 'bmrb_id':
                kk = '{}-{}-{}-{}-{}-{}'.format(amide_info[0], amide_info[1], amide_info[2], amide_info[3], aro1_info[0],
                                                aro1_info[1])
                #if kk in k.keys():# and abs(float(amide_info[5]))>2 and k[kk]!='unlikely':
                if int(amide_info[7])==1 and -15.0 < float(amide_info[5]) < 15.0:# and float(aro1_info[2])<3.0:
                        #(and float(amide_info[5])>sigma_cutoff or float(amide_info[5]) < -sigma_cutoff) \
                        #and float(aro1_info[2])<8.0:
                    l='{}/{}/{}/{}/{}/{}'.format(amide_info[0],amide_info[1],amide_info[2],amide_info[3],aro1_info[0],aro1_info[1])
                    pdb_subset.append(amide_info[0])
                    if bid!=amide_info[1]:
                        sw=get_software_info(amide_info[1])
                        bid = amide_info[1]
                    if len(sw)>0:
                        if True in [True for i in sw if 'rosetta' in i.lower()]:
                            sw_info.append('CS_ROSETTA')
                        elif True in [True for i in sw if 'gamdy' in i.lower()]:
                            sw_info.append('CS-GAMDY')
                        elif True in [True for i in sw if 'plor' in i.lower()]:
                            sw_info.append('X-PLOR')
                        elif True in [True for i in sw if 'cyana' in i.lower()]:
                            sw_info.append('CYANA')
                        elif True in [True for i in sw if 'aria' in i.lower()]:
                            sw_info.append('ARIA')
                        elif True in [True for i in sw if 'felix' in i.lower()]:
                            sw_info.append('FELIX')
                        elif True in [True for i in sw if 'amber' in i.lower()]:
                            sw_info.append('AMBER')
                        elif True in [True for i in sw if 'talos' in i.lower()]:
                            sw_info.append('TALOS')
                        elif True in [True for i in sw if 'nmrfam-sparky' in i.lower()]:
                            sw_info.append('NMRFAM-SPARKY')
                        elif True in [True for i in sw if 'cns' in i.lower()]:
                            sw_info.append('CNS')
                        else:
                            sw_info.append('Other')
                            #print (",".join(sw))
                    else:
                        sw_info.append('No info')
                    print (amide_info[1],sw_info[-1])
                    lab.append(l)
                    mean_d.append(float(aro1_info[2]))
                    std_d.append(float(aro1_info[3]))
                    z.append(float(amide_info[5]))
                    cs.append(float(amide_info[4]))
                    if float(amide_info[5])< -sigma_cutoff:
                        sigma.append('z < -{}'.format(sigma_cutoff))
                    elif float(amide_info[5]) > sigma_cutoff:
                        sigma.append('z > {}'.format(sigma_cutoff))
                    else:
                        sigma.append('-{} < z < {}'.format(sigma_cutoff,sigma_cutoff))
                    if float(aro1_info[2])<8.0:
                        d_cutoff.append('Distance<8A')
                    else:
                        d_cutoff.append('Full data')
                    # zrange=numpy.arange(-15.0,15.0,1.0)
                    # for i in range(len(zrange)):
                    #     if i<len(zrange)-1:
                    #         if zrange[i+1] > float(amide_info[5]) >= zrange[i]:
                    #             ccc.append('{}sigma>z>={}sigma'.format(zrange[i+1],zrange[i]))
            # fig = px.scatter(x=mean_d,y=z,color=color,opacity=0.5,
            #                  #facet_col=color,
            #                  labels=dict(x='Mean distance', y='Z-Score'))
            # #fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
            # fig.update_layout(legend=dict(title=''))
        x1=[]
        y1=[]
        z1=[]
        for i in range(len(sw_info)):
            if sw_info[i] == 'CS_ROSETTA':
                x1.append(mean_d[i])
                y1.append(z[i])
                z1.append(sw_info[i])

        fig = px.scatter(x=x1,y=y1,
                         #color=sw_info,
                         #hover_name=lab,
                           # opacity=0.6,
                           #facet_row=d_cutoff,
                           #facet_col=color,
                           # histnorm='probability',
                           labels=dict(x='Distance (Å)', y='Z-score'),
                            range_y=[-5,5],
                            range_x=[0,100],
                           #category_orders={"facet_col": ["TYR", "PHE", "TRP", "HIS"]}
                           )
        #fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_traces(marker=dict(size=3))
        #fig.update_xaxes(matches=None)
        #fig.update_xaxes(showticklabels=True, row=1)
        #fig.update_xaxes(showticklabels=True, row=2)
        #fig.update_layout(showlegend=False)
        fig.update_layout(showlegend=False, font=(dict(family='Arial', size=15, color='black')))
        outfile = 'dist_vs_z_sw_3'
        fig.write_html('/Users/kumaran/Documents/rcs_manuscript/{}.html'.format(outfile))
        fig.write_image('/Users/kumaran/Documents/rcs_manuscript/{}.jpeg'.format(outfile), width=1200, height=800)
        fig.write_image('/Users/kumaran/Documents/rcs_manuscript/{}.pdf'.format(outfile), width=1200, height=800)
        fig.show()
        print (len(set(pdb_total)))
        print (len(set(pdb_subset)))


if __name__=="__main__":
    csv_file = sys.argv[1]
    plot_software_info(csv_file)