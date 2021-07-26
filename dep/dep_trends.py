import json
import requests
import csv
import plotly.graph_objects as go

url = 'https://api.bmrb.io/v2/meta/release_statistics'
r = requests.get(url).json()
release_dict = r['release_information']
years = []
total = []

for year in release_dict:
    release = release_dict[year]['structure_release_in_year'] # Total structural depositions
    if release['total'] != 0:
        years.append(year)
        total.append(release['total'])

csr_filename = 'CS-Rosetta_Entries.csv'
csr_by_year = {}
csr_list = [0 for year in years]
with open(csr_filename, 'r') as csfile:
    csv_reader = csv.reader(csfile, delimiter=',')
    next(csv_reader)
    for row in csv_reader:
        year = row[-1]
        for i in range(len(years)):
            if years[i] == year:
                year_index = i
        csr_list[year_index] += 1

csr_sub = [
    0, 0, 0, 0, 0, 0, 0, 9, 621, 571, 676, 1195, 975, 647,
    1098, 931, 908, 658, 332
] # taken from https://csrosetta.bmrb.io/ ; the final number will change as 2021 progresses
csr_totals = []

fig = go.Figure()
fig.add_trace(go.Scatter(
    x=years[:-2], y=total[:-2], mode='lines+markers', name='Total Structure Depositions'
))
fig.add_trace(go.Scatter(
    x=years[:-2], y=csr_list[:-2], mode='lines+markers', name='CS-Rosetta Depositions'
))
fig.add_trace(go.Scatter(
    x=years[:-2], y=csr_sub[:-2], mode='lines+markers', name='CS-Rosetta Runs'
))

fig.update_layout(
    xaxis_title='Year',
    yaxis_title='Instances in Year',
    font = dict(family='Arial', size=18),
    legend=dict(y=0.95, x=0.05)
)

#fig.show(renderer="firefox")
fig.write_image("../images/dep_plot.pdf")



