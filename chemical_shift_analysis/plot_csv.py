import csv
import plotly.express as px

def plot_d_vs_z(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z=[]
        d=[]
        ang=[]
        sang=[]
        lab=[]
        ar_res=[]
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0],row[1],row[4],row[6],row[9],row[11],row[7],row[8],row[12])
            if abs(float(row[8])) < 20 and int(row[2])==1 and int(row[3])==1:
                d.append(float(row[12]))
                z.append(float(row[8]))
                ang.append(float(row[17]))
                sang.append(float(row[22]))
                ar_res.append(row[11])
                lab.append(info)

    fig = px.scatter(x=d,y=z,color=ar_res,hover_name=lab,
                     labels=dict(x='Distance (Å)', y='Z-score'),)
    fig.show()


def plot_d_vs_solidangle(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z=[]
        d=[]
        ang=[]
        sang=[]
        lab=[]
        ar_res=[]
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0],row[1],row[4],row[6],row[9],row[11],row[7],row[8],row[12])
            if abs(float(row[8])) < 20 and int(row[2]) == 1 and int(row[3]) == 1:
                d.append(float(row[12]))
                z.append(float(row[8]))
                ang.append(float(row[17]))
                sang.append(float(row[22]))
                ar_res.append(row[11])
                lab.append(info)
    fig = px.scatter(x=d,y=sang,color=ar_res,hover_name=lab,
                     labels=dict(x='Distance (Å)', y='Solid angle'),)
    fig.show()

def plot_azimuthal_vs_z(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z = []
        d = []
        ang = []
        sang = []
        lab = []
        ar_res = []
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0], row[1], row[4], row[6], row[9], row[11], row[7],
                                                       row[8], row[12])
            if abs(float(row[8])) < 20 and int(row[2])==1 and int(row[3])==1:
                d.append(float(row[12]))
                z.append(float(row[8]))
                ang.append(float(row[17]))
                sang.append(float(row[22]))
                ar_res.append(row[11])
                lab.append(info)
    fig = px.scatter(x=ang, y=z, color=ar_res, hover_name=lab,
                     labels=dict(x='Azimuthal angle', y='Z-score'), )
    fig.show()

def plot_solidangle_vs_z(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z = []
        d = []
        ang = []
        sang = []
        lab = []
        ar_res = []
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0], row[1], row[4], row[6], row[9], row[11], row[7],
                                                       row[8], row[12])
            if abs(float(row[8])) < 20 and int(row[2])==1 and int(row[3])==1:
                d.append(float(row[12]))
                z.append(float(row[8]))
                ang.append(float(row[17]))
                sang.append(float(row[22]))
                ar_res.append(row[11])
                lab.append(info)
    fig = px.scatter(x=sang, y=z, color=ar_res, hover_name=lab,
                     labels=dict(x='Solid angle', y='Z-score'), )
    fig.show()

def plot_3d(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z = []
        d = []
        ang = []
        sang = []
        lab = []
        ar_res = []
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0], row[1], row[4], row[6], row[9], row[11], row[7],
                                                       row[8], row[12])
            if abs(float(row[8])) < 20 and int(row[2])==1 and int(row[3])==1:
                d.append(float(row[12]))
                z.append(float(row[8]))
                ang.append(float(row[17]))
                sang.append(float(row[22]))
                ar_res.append(row[11])
                lab.append(info)
    fig = px.scatter_3d(x=d, y=ang, z=z, color=ar_res, hover_name=lab,
                     labels=dict(x='Distance (Å)', z='Z-score', y='Azimuthal angle'), )
    fig.show()

def plot_azimithal_solid(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z = []
        d = []
        ang = []
        sang = []
        lab = []
        ar_res = []
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0], row[1], row[4], row[6], row[9], row[11], row[7],
                                                       row[8], row[12])
            if abs(float(row[8])) < 20 and int(row[2])==1 and int(row[3])==1:
                d.append(float(row[12]))
                z.append(float(row[8]))
                ang.append(float(row[17]))
                sang.append(float(row[22]))
                ar_res.append(row[11])
                lab.append(info)
    fig = px.scatter(x=ang, y=sang,  color=ar_res, hover_name=lab,
                     labels=dict( x='Azimuthal angle',y='Solid angle') )
    fig.show()


def azimuthal_hist(csv_file):
    with open(csv_file, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        z = []
        d = []
        ang = []
        sang = []
        lab = []
        ar_res = []
        sigma = []
        z2 = []
        sigma_cutoff = 3.0
        for row in csvreader:
            info = '{}-{}-{}-{}-{}-{}/{}/{}/{}'.format(row[0], row[1], row[4], row[6], row[9], row[11], row[7], row[8],
                                                       row[12])
            if abs(float(row[8])) < 20 and int(row[2]) == 1 and int(row[3]) == 1 and abs(float(row[8]))>sigma_cutoff:
                d.append(float(row[12]))
                z.append(float(row[8]))
                ang.append(float(row[17]))
                sang.append(float(row[22]))
                ar_res.append(row[11])
                lab.append(info)
                if float(row[8]) < -sigma_cutoff:
                    sigma.append('Upfield (Z < -{})'.format(sigma_cutoff))
                elif float(row[8]) > sigma_cutoff:
                    sigma.append('Downfield (Z > {})'.format(sigma_cutoff))
                else:
                    sigma.append('-{}<= Z-score <= {}'.format(sigma_cutoff, sigma_cutoff))

                zrange = [-20.0,-15.0, -5.0, -4.0, -3.0, -2.0, 2.0, 3.0, 4.0, 5.0, 15.0,20.0]
                for i in range(len(zrange)):
                    if i < len(zrange) - 1:
                        if zrange[i + 1] > float(row[8]) >= zrange[i]:
                            z2.append('${}\sigma \gt Z \ge {} \sigma$'.format(zrange[i + 1], zrange[i]))
    bin_width = 3
    bins = int(round((max(sang) - min(sang)) / bin_width, 0))
    fig = px.histogram(x=sang, color=z2, facet_row=sigma,
                       facet_col=ar_res,
                       nbins=bins,
                       labels={"x": 'Azimuth angle (°)', "count": 'Number'},
                       category_orders={"facet_col": ["TYR", "PHE", "TRP", "HIS"]}
                       )#.update_layout(yaxis_title='Count', yaxis5_title='Count')

    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    # fig.update_layout(barmode='overlay')
    fig.update_layout(
        legend=dict(
            yanchor="top",
            y=0.99,
            title='',
            xanchor="left",
            x=0.001, bgcolor='rgba(0,0,0,0)'
        ), font=(dict(family='Arial', size=15, color='black')))
    fig.show()
if __name__ == "__main__":
    # plot_d_vs_z('/Users/kumaran/Downloads/output/output.txt')
    # plot_azimuthal_vs_z('/Users/kumaran/Downloads/output/output.txt')
    # plot_solidangle_vs_z('/Users/kumaran/Downloads/output/output.txt')
    # plot_azimithal_solid('/Users/kumaran/Downloads/output/output.txt')
    azimuthal_hist('/Users/kumaran/Downloads/output/output.txt')