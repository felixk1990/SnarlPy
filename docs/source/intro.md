About goflow: 
A repository holding the structure for the 'intertwined' package as well
examples and galleries.

##  Introduction
The module 'intertwined' is a python packages encompassing a set of class and
method implementations for networkx datatypes, in order to calculate linking
numbers of spatially intertwined networks and near optimal cuts to
topologically unlink them. Used and explained in the publication: ENTER_REF
<br>

##  Installation
Via PyPi
```
pip install intertwined
```
Via Cloning
```
git clone https://git.mpi-cbg.de/kramer/entanglement-analysis.git .
pip install ./entanglement-analysis/
```
##  Usage
For theory and algorithm details see publication: ENTER_REF
Let's create an intertwined, spatial network with the package's internal
generator (here two ladders) and calculate the linking number matrices 'lk_mat'
(cycle space) and 'p'(edge space). When we have done so compute the the priority
matrices and plot the matrices diagnonal values onto the respective graphs'
edges (visualising the important edges for intertwinedness).
```
import numpy as np
import intertwined.edgePriority as itwe
import intertwined.tangledGenerators as tg

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

num_periods = 2
D = tg.createLabelCatenation(num_periods)
graph_sets = [k.G for k in D.layer]
p, lk_mat = itwe.getEdgeLinkageOperator(graph_sets)
fig = update_priorityPlot(D, p)
fig.show()
```
![network](https://git.mpi-cbg.de/kramer/entanglement-analysis/-/raw/main/gallery/main/dualLadderShift_0.png)
![prio1](https://git.mpi-cbg.de/kramer/entanglement-analysis/-/raw/main/gallery/main/lambdaSQ_10.png)
![prio2](https://git.mpi-cbg.de/kramer/entanglement-analysis/-/raw/main/gallery/main/lambdaSQ_20.png)

The package also allows you to directly compute cut sets, to find the best way
to topologically disentangle the networks. Calling this bit of code will
repeatedly compute the priority matrices for both networks and remove edges of
highest priority one by one (as displayed on the carton below).
```
import intertwined.sampling as itws
import intertwined.edgePriority as itwe
import tangledGenerators as tg

num_periods = 2
D = tg.createLabelHexagonHopfed_V4()
graph_sets = [k.G for k in D.layer]
init_cut_sets = [graph_sets[:], graph_sets[::-1]]
cut_lists = []
for ics in init_cut_sets:
    cut_lists.append(itwe.cuttingEdgeAlgorithm(*ics))
```
![cuts](https://git.mpi-cbg.de/kramer/entanglement-analysis/-/raw/main/gallery/main/cuttingEdgeAlgorithm.png)
##  Requirements
```
networkx==2.5
numpy==1.21.0
scipy==1.7.3
kirchhoff==0.2.7
```

## Acknowledgement
```intertwined``` written by Felix Kramer
