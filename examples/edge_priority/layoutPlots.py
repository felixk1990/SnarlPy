import networkx as nx
import numpy as np
import scipy
import os.path as op

import plotly.graph_objects as go
import iop.pgf_figure as pf
import iop.init_IO as initIO

pathOutput = './plots_cuttingEdge'

def plotLambdaMaps(p, set_number):

    myMaps=['BuGn','RdPu']
    opt =dict(width="140%",
            height="5%",
            loc='lower center',
            borderpad=-4
         )

    I = [np.dot(p,p.T), np.dot(p.T,p)]
    for i in range(2):

        fig, ax = pf.newfig(0.4)
        tag = r'$\Vert S_{'+str(i+1)+',ij} \Vert$'

        im = ax.imshow(np.absolute(I[i]),cmap=myMaps[i])
#         axins = inset_axes(ax,**opt)
#         cbar = fig.colorbar(im, cax=axins , orientation="horizontal")
        cbar = fig.colorbar(im, ax=ax )
        ax.set_title(tag)

        axtag =r'$E_'+str(i+1)+'$'
        ax.set_xlabel(axtag)
        ax.set_ylabel(axtag)

#         ticks = [*range(len(I[i][0]))][::4]
#         ax.set_xticks(ticks)
#         ax.set_yticks(ticks)

        file = 'pyplotLambdaSQ_'+str(i+1)+str(set_number)
        filename=op.join(pathOutput, file)
        pf.savefig(filename)

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

def export_priorityPlot(D, p, file):

    global pathOutput

    fig = update_priorityPlot(D, p)

    filename=op.join(pathOutput, file)
    fig.write_html(filename +".html" )
    fig.write_image(filename +".pdf" )
