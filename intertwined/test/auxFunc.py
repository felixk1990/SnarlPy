# @Author: Felix Kramer <felix>
# @Date:   2022-07-28T10:36:11+02:00
# @Email:  felixuwekramer@proton.me
# @Filename: auxFunc.py
# @Last modified by:   felix
# @Last modified time: 2022-07-28T10:38:40+02:00
import networkx as nx
import numpy as np
import kirchhoff.circuit_init as kc
import kirchhoff.circuit_dual as kd


def contruct_link_0():

    D = kc.initialize_circuit_from_crystal('square', (1, 2))

    # unlink
    phi = np.linspace(0, 1, num=50) * 2. * np.pi
    R = 0.5
    G0 = nx.Graph()
    for i, p in enumerate(phi[:-1]):
        XYZ = (3., 2.+np.cos(p), np.sin(p))
        G0.add_node(i, pos=np.array(XYZ)*R)

    for i, p in enumerate(phi[:-2]):
        G0.add_edge(i, i+1)

    G0.add_edge(len(phi)-2, 0)
    K = kd.DualCircuit([kc.Circuit(G0), D])

    return K


def contruct_link_1():

    D = kc.initialize_circuit_from_crystal('square', (1, 2))

    # simple link
    phi = np.linspace(0, 1, num=50) * 2. * np.pi
    R = 0.5
    G1 = nx.Graph()
    for i, p in enumerate(phi[:-1]):
        XYZ = (np.cos(p), 1., np.sin(p))
        G1.add_node(i, pos=np.array(XYZ)*R)

    for i, p in enumerate(phi[:-2]):
        G1.add_edge(i, i+1)

    G1.add_edge(len(phi)-2, 0)
    K = kd.DualCircuit([kc.Circuit(G1), D])

    return K


def contruct_link_2():

    D = kc.initialize_circuit_from_crystal('square', (1, 2))

    # double link
    phi = np.linspace(0, 1, num=50) * 2. * np.pi
    R = 0.5
    G2 = nx.Graph()
    for i, p in enumerate(phi[:-1]):
        XYZ = (np.cos(2.*p+np.pi), 1+p/(np.pi), np.sin(2.*p+np.pi))
        G2.add_node(i, pos=np.array(XYZ)*R)

    for i, p in enumerate(phi[:-2]):
        G2.add_edge(i, i+1)

    G2.add_node(len(phi)-1, pos=np.array((-1, 1, 0)))
    G2.add_edge(len(phi)-2, len(phi)-1)
    G2.add_edge(0, len(phi)-1)

    K = kd.DualCircuit([kc.Circuit(G2), D])

    return K


def contruct_link_3():

    D = kc.initialize_circuit_from_crystal('square', (1, 2))

    # simple link no 2
    phi = np.linspace(0, 1, num=50) * 2. * np.pi
    R = 0.5
    G0 = nx.Graph()
    for i, p in enumerate(phi[:-1]):
        XYZ = (1., 2.+np.cos(p), np.sin(p))
        G0.add_node(i, pos=np.array(XYZ)*R)

    for i, p in enumerate(phi[:-2]):
        G0.add_edge(i, i+1)

    G0.add_edge(len(phi)-2, 0)

    K = kd.DualCircuit([kc.Circuit(G0), D])

    return K


constr_gen = [
    contruct_link_0,
    contruct_link_1,
    contruct_link_2,
    contruct_link_3,
]
