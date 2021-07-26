import plotly.graph_objects as go
import pandas as pd

separate_energies = 14.6385 + 1.75457 # the VdW self interaction for PHE and ALA, respectively
# this can be found by running energy calculations in mosart
# of isolated PHE and isolated ALA

df = pd.read_csv('energies_fine.csv') # the file where the assembly VdW energies have been stored
# according to the distance from the ALA amide H to the PHE ring center
df['VdW'] = df['VdW'] - separate_energies # subtracting the self interaction


#plotting
fig = go.Figure()
fig.add_trace(go.Scatter(
    x=df['Distance'],
    y=df['VdW'],
    mode='markers+lines',
    marker=dict(
        color='cadetblue',
        size=16,
        line=dict(width=2)
    ),
    textfont_size=18,
    textfont_family="Arial"
))
fig.add_hline(y=3.3) # the magnitude of the attractive potential given by Levitt and Perutz 1988
fig.update_layout(
    xaxis_title = "Distance Between Ring Center and Amide Nitrogen",
    yaxis_title = "VdW Energy (kcal/mol)",
    title_x=0.5,
    font=dict(family="Arial",size=18)
)
fig.update_yaxes(range=[-5, 30])

fig.show(renderer='firefox')
fig.write_image("plot_vdw.pdf")
