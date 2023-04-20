import networkx as nx
import numpy as np
import scipy
import os.path as op

import plotly.graph_objects as go

def plotLambdaMaps(p):

    myMaps=['BuGn','RdPu']
    opt =dict(width="140%",
            height="5%",
            loc='lower center',
            borderpad=-4
         )
    cmap = ['BuGn', 'RdPu']
    I = [np.dot(p,p.T), np.dot(p.T,p)]
    for i in range(2):
        
        kwargs = dict(
            autosize = False,
            height = 400,
            width = 400,
            coloraxis = {'colorscale':cmap[i]},
            title={
                'text':  r'$\Vert S_{'+str(i+1)+',ij} \Vert$',
                'y':.9,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'},
            xaxis_title=r'$E_'+str(i+1)+'$',
            yaxis_title=r'$E_'+str(i+1)+'$',
        )
        
        fig = go.Figure()    
        fig.add_trace(go.Heatmap(z=np.absolute(I[i]), coloraxis="coloraxis"))
        fig.update_layout(**kwargs)
        fig.show()

def update_plot(D):

    kwargs = {
        'axis': False
    }

    fig = D.plot_circuit('centrality', **kwargs)
    fig.show()

    return fig

def update_priorityPlot(D, p):

    colors = []
    priority_mat = [np.dot(p,p.T), np.dot(p.T,p)]

    for i, k in enumerate(D.layer):
        k.edges['priority'] = np.round( np.diagonal(priority_mat[i]),10)
        colors.append(k.edges['priority'])
    kwargs = {
        'color_nodes': ['#030512','#030512'],
        'color_edges': colors,
        'colormap': ['BuGn','RdPu'],
        'axis': True
    }

    fig = D.plot_circuit('priority', **kwargs)
    fig.show()

    return fig
