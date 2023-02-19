# @Author:  Felix Kramer
# @Date:   2021-11-27T20:10:14+01:00
# @Email:  kramer@mpi-cbg.de
# @Project: go-with-the-flow
# @Last modified by:    Felix Kramer
# @Last modified time: 2021-11-28T21:53:35+01:00
# @License: MIT



import sys
import os.path as op
import networkx as nx
import numpy as np
import scipy

import kirchhoff.circuit_init as kfi
import kirchhoff.circuit_dual as kcd
import intertwined.sampling as itws
import intertwined.edgePriority as itwe

def createLabelCatenation(num_periods):

    nx_type = 'catenation'
    D = kcd.initialize_dual_from_catenation(dual_type=nx_type, num_periods=num_periods)

    for i, k in enumerate(D.layer):

        k.nodes['label'] = [n for n in k.G.nodes()]
        k.edges['label'] = [e for i, e in enumerate(k.G.edges())]

    return D

def createLabelLaves(num_periods):

    nx_type = 'laves'
    D = kcd.initialize_dual_circuit_from_minsurf(dual_type=nx_type, num_periods=num_periods)

    for i, k in enumerate(D.layer):

        k.nodes['label'] = [n for n in k.G.nodes()]
        k.edges['label'] = [e for i, e in enumerate(k.G.edges())]

    return D

def createLabelLavesPruned(num_periods, pruning, idx):

    nx_type = 'laves'
    D1 = kcd.initialize_dual_circuit_from_minsurf(dual_type=nx_type, num_periods=num_periods)

    for i, k in enumerate(D1.layer):

        if i == idx:
            for p in pruning:
                k.G.remove_edge(*p)

            list_n = [n for n in k.G.nodes() if k.G.degree(n)==0]
            for n in list_n:
                k.G.remove_node(n)

    D2=kcd.initialize_dual_circuit_from_networkx(D1.layer[0].G,D1.layer[1].G,[])

    for i, k in enumerate(D2.layer):
        k.nodes['label'] = [n for n in k.G.nodes()]
        k.edges['label'] = [e for i, e in enumerate(k.G.edges())]

    return D2

def shiftLadder(D, step):

    for i, k in enumerate(D.layer):
        if i == 0:

            for n in k.G.nodes():

                k.G.nodes[n]['pos'] = np.add(k.G.nodes[n]['pos'], step)

bouncerList = [(53,42), (54,55), (46,57), (44,55),(40,51), (38,49), (36,47), (47,48), (48,49), (49,50), (50,51), (51,52), (52,53), (53,54), (55,56), (56,57) , (45,46)]

def createHexagonal(num_periods=2):

    K=kfi.initialize_circuit_from_crystal('hexagonal',num_periods)

    for n in K.G.nodes():
        p=K.G.nodes[n]['pos']
        K.G.nodes[n]['pos']=np.array([*p,0.])

    return K

def createSquare(num_periods=2):

    K=kfi.initialize_circuit_from_crystal('square',num_periods)
    #
    # for n in K.G.nodes():
    #     p=K.G.nodes[n]['pos']
    #     K.G.nodes[n]['pos']=np.array([*p,0.])

    return K

def createHexagon(num_periods=2):

    K=kfi.initialize_circuit_from_crystal('hexagonal',num_periods)

    for n in K.G.nodes():
        p=K.G.nodes[n]['pos']
        K.G.nodes[n]['pos']=np.array([*p,0.])

    global bouncerList
    for bl in bouncerList:
        K.G.remove_edge(*bl)

    list_n = [n for n in K.G.nodes() if K.G.degree(n)==0]
    for n in list_n:
        K.G.remove_node(n)

    return K

def createLabelHexagonHopfed_V1():

    K = createHexagon()
    ##### generate intertwined curves
    phi=np.linspace(0,1,num=20)*2.*np.pi
    R=5

    # unlink
    G0=nx.Graph()
    for i,p in enumerate(phi[:-1]):
        XYZ=(2.5,.5+R*np.cos(p),np.sin(p)*R)
        G0.add_node(i,pos=np.array(XYZ))

    for i,p in enumerate(phi[:-2]):
        G0.add_edge(i,i+1)

    G0.add_edge(len(phi)-2,0)

    D=kcd.initialize_dual_circuit_from_networkx(G0,K.G,[])

    return D

def createLabelHexagonHopfed_V2():

    K = createHexagon()
    ##### generate intertwined curves
    phi=np.linspace(0,1,num=20)*2.*np.pi
    R=2

    # unlink
    G0=nx.Graph()
    for i,p in enumerate(phi[:-1]):
        XYZ=(2.5,3.5+R*np.cos(p),np.sin(p)*R)
        G0.add_node(i,pos=np.array(XYZ))

    for i,p in enumerate(phi[:-2]):
        G0.add_edge(i,i+1)

    G0.add_edge(len(phi)-2,0)

    D=kcd.initialize_dual_circuit_from_networkx(G0,K.G,[])

    return D

def createLabelHexagonHopfed_V3():

    K = createHexagon()
    ##### generate intertwined curves
    phi=np.linspace(0,1,num=20)*2.*np.pi
    R=3.5

    # unlink
    G0=nx.Graph()
    for i,p in enumerate(phi[:-1]):
        XYZ=(2.5,5.+R*np.cos(p),np.sin(p)*R)
        G0.add_node(i,pos=np.array(XYZ))

    for i,p in enumerate(phi[:-2]):
        G0.add_edge(i,i+1)

    G0.add_edge(len(phi)-2,0)

    D=kcd.initialize_dual_circuit_from_networkx(G0,K.G,[])

    return D

def createLabelHexagonHopfed_V4():

    K = createHexagon()
    ##### generate intertwined curves
    phi=np.linspace(0,1,num=20)*2.*np.pi
    R=2.5

    # unlink
    G0=nx.Graph()
    for i,p in enumerate(phi[:-1]):
        XYZ=(2.5,4.25+R*np.cos(p),np.sin(p)*R)
        G0.add_node(i,pos=np.array(XYZ))

    for i,p in enumerate(phi[:-2]):
        G0.add_edge(i,i+1)

    G0.add_edge(len(phi)-2,0)

    D=kcd.initialize_dual_circuit_from_networkx(G0,K.G,[])

    return D

def createLabelHexagonHopfed_V5():

    K = createHexagon(4)
    ##### generate intertwined curves
    phi=np.linspace(0,1,num=20)*2.*np.pi
    R=2.5

    # unlink
    G0=nx.Graph()
    for i,p in enumerate(phi[:-1]):
        XYZ=(5.5,4.25+R*np.cos(p),np.sin(p)*R)
        #5.75
        G0.add_node(i,pos=np.array(XYZ))

    for i,p in enumerate(phi[:-2]):
        G0.add_edge(i,i+1)

    G0.add_edge(len(phi)-2,0)

    D=kcd.initialize_dual_circuit_from_networkx(G0,K.G,[])

    return D

def createLabelDiamond(num_periods):

    nx_type = 'diamond'
    D = kcd.initialize_dual_circuit_from_minsurf(dual_type=nx_type, num_periods=num_periods)

    for i, k in enumerate(D.layer):

        k.nodes['label'] = [n for n in k.G.nodes()]
        k.edges['label'] = [e for i, e in enumerate(k.G.edges())]

    return D

def createLabelSimple(num_periods):

    nx_type = 'simple'
    D1 = kcd.initialize_dual_circuit_from_minsurf(dual_type=nx_type, num_periods=num_periods)
    for dl in D1 .layer:
        list_n=list(dl.G.nodes())
        idx =[n for n in list_n if dl.G.degree(n)==3]
        for i in idx:
            dl.G.remove_node(i)

    D2=kcd.initialize_dual_circuit_from_networkx(D1.layer[0].G,D1.layer[1].G,[])

    for i, k in enumerate(D2.layer):

        k.nodes['label'] = [n for n in k.G.nodes()]
        k.edges['label'] = [e for i, e in enumerate(k.G.edges())]

    return D2

def createLabelSimplePruned(num_periods, pruning, idx):

    nx_type = 'simple'
    D1 = kcd.initialize_dual_circuit_from_minsurf(dual_type=nx_type, num_periods=num_periods)

    for dl in D1 .layer:
        list_n=list(dl.G.nodes())
        idx =[n for n in list_n if dl.G.degree(n)==3]
        for i in idx:
            dl.G.remove_node(i)

    for i, k in enumerate(D1.layer):

        if i == idx:
            for p in pruning:
                k.G.remove_edge(*p)

            list_n = [n for n in k.G.nodes() if k.G.degree(n)==0]
            for n in list_n:
                k.G.remove_node(n)

    D2=kcd.initialize_dual_circuit_from_networkx(D1.layer[0].G,D1.layer[1].G,[])

    for i, k in enumerate(D2.layer):
        k.nodes['label'] = [n for n in k.G.nodes()]
        k.edges['label'] = [e for i, e in enumerate(k.G.edges())]

    return D2

def createLabelHexagonStack(num_periods):

    K1 = kfi.initialize_circuit_from_crystal('hexagonal',num_periods+1)
    K2a = kfi.initialize_circuit_from_crystal('hexagonal',num_periods)
    K2b = kfi.initialize_circuit_from_crystal('hexagonal',num_periods)

    push = []
    for i, k in enumerate([K2a,K1,K2b]):

        k.nodes['label'] = [n for n in k.G.nodes()]
        k.edges['label'] = [i for i, e in enumerate(k.G.edges())]
        for n in k.G.nodes():
#             k.G.nodes[n]['pos'] = np.add(k.G.nodes[n]['pos'], [0,0,i])
            k.G.nodes[n]['pos'] = [*k.G.nodes[n]['pos'], i]
            if i == 1:
                p = k.G.nodes[n]['pos']
                shift = np.array([-1, -1*np.sqrt(3)/2., 0])*2.
                q = np.add(shift,p)
                k.G.nodes[n]['pos']=q

                push.append(q[:-1])

    H = nx.Graph()
    for i, g in enumerate([K2a.G,K2b.G]):

        for n in g.nodes():

            p = g.nodes[n]['pos']
            H.add_node(tuple(p), pos=p)

            if i > 0:
                q = np.subtract(p, [0, 0, 2])
                d = np.linalg.norm(np.subtract(q[:-1], push), axis=1)
#                 print(d)
                if np.all(d > 0.001):
                    H.add_edge(tuple(p), tuple(q))

        for e in g.edges():

            p = g.nodes[e[0]]['pos']
            q = g.nodes[e[1]]['pos']
            H.add_edge(tuple(p), tuple(q))

    D = kcd.initialize_dual_circuit_from_networkx(K1.G, H,[])
    for l in D.layer:
         l.nodes['pos'] = [l.G.nodes[n]['pos'] for n in l.G.nodes()]

    return D

def createRandomHopfed():

    # K = kf.initialize_circuit_from_random('voronoi_planar', 30, 1)
    # nx.write_gpickle(K.G,'./pickle_rick/nx_random.pkl')
    G=nx.read_gpickle('./pickle_rick/nx_random.pkl')
    for n in G.nodes():
        p=G.nodes[n]['pos']
        G.nodes[n]['pos']=np.array([*p,0.])

    ##### generate intertwined curves
    phi=np.linspace(0,1,num=20)*2.*np.pi
    R=0.3

    # unlink
    G0=nx.Graph()
    for i,p in enumerate(phi[:-1]):
        XYZ=(0.5,0.5+R*np.cos(p),np.sin(p)*R)
        G0.add_node(i,pos=np.array(XYZ))
    for i,p in enumerate(phi[:-2]):
        G0.add_edge(i,i+1)

    G0.add_edge(len(phi)-2,0)

    D=kcd.initialize_dual_circuit_from_networkx(G0,G,[])

    return D

def createRandomMeshed():

    # K = kfi.initialize_circuit_from_random('voronoi_planar', 7, 5)
    # nx.write_gpickle(K.G,'./pickle_rick/nx_random.pkl')
    G1 = nx.read_gpickle('./pickle_rick/nx_random.pkl')
    for n in G1.nodes():
        p = G1.nodes[n]['pos']
        G1.nodes[n]['pos'] = np.array([*p,0.])

    # G2 = nx.read_gpickle('./pickle_rick/nx_random.pkl')

    # rx = 1.
    # ry =rx
    # dz = -1.
    # for n in G2.nodes():
    #     p = np.add(G2.nodes[n]['pos'], [rx*np.random.rand(),ry*np.random.rand()])
    #     G2.nodes[n]['pos'] = np.array([*p, dz*(1-2.*np.random.rand())])
    # nx.write_gpickle(G2,'./pickle_rick/nx_random2.pkl')
    G2 = nx.read_gpickle('./pickle_rick/nx_random2.pkl')

    # D = kcd.initialize_dual_circuit_from_networkx(G1,G2,[])
    import kirchhoff.circuit_dual as kcd
    import kirchhoff.init_dual as kci
    d = kci.NetworkxDual()
    d.layer = [G1, G2]
    D = kcd.DualCircuit(kcd.construct_from_graphSet(d))

    return D

def createRandomMeshed3D():

    K = kfi.initialize_circuit_from_random('voronoi_volume', 5, 5)
    nx.write_gpickle(K.G,'./pickle_rick/nx_random3D1.pkl')

    G1 = nx.read_gpickle('./pickle_rick/nx_random3D1.pkl')
    G2 = nx.read_gpickle('./pickle_rick/nx_random3D1.pkl')

    scale = 1.
    for n in G2.nodes():
        sign = (-1)**(int(np.random.rand()))
        p = np.add(G2.nodes[n]['pos'], scale*np.random.rand(3)*sign)
        G2.nodes[n]['pos'] = p

    nx.write_gpickle(G2, './pickle_rick/nx_random3D2.pkl')
    # D = kcd.initialize_dual_circuit_from_networkx(G1,G2,[])
    import kirchhoff.circuit_dual as kcd
    import kirchhoff.init_dual as kci
    d = kci.NetworkxDual()
    d.layer = [G1, G2]
    D = kcd.DualCircuit(kcd.construct_from_graphSet(d))

    return D
