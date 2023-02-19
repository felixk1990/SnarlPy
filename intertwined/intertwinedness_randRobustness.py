# @Author: Felix Kramer <felix>
# @Date:   2022-07-27T15:50:36+02:00
# @Email:  felixuwekramer@proton.me
# @Filename: intertwinedness_randRobustness.py
# @Last modified by:   felix
# @Last modified time: 2022-07-28T18:42:59+02:00
import networkx as nx
import numpy as np
import random as rd
import intertwined.sampling as itws


def sample_linkage_robustness_onesided(ref_graph, cut_graph, repNum=1):
    """
    sampling through the cycle space of a seeminly intergraph system.
    Testing whether sampled cycles are intertwined with each other.

    Args:
        ref_graph (list):\n
            The reference graph, which is setting the frame to be disentangled
            from.
        cut_graph (list):\n
            The graph for which the optimal cut set if searched.

    Returns:
        ndarray:\n
            The relative cutting numbers for each sampling.
        list:\n
            The list of optimal cut set for the cut_graph.

    """
    robustnessSets = []
    cutSets = []
    nullitySets = []

    for i in range(repNum):

        # generate basis vectors and corresponding graph matrices
        graph_sets = [cut_graph, ref_graph]
        cyc_nx_base = [itws.calc_cycle_basis(G) for G in graph_sets]
        nullity = np.array([itws.calc_nullity(G)-1. for G in graph_sets])
        nullitySets.append(nullity[0])

        # compute linkage of basis vectors
        numeric_res = itws.basis_linkage(*cyc_nx_base)
        res = np.round(
            np.fromiter(numeric_res[(0, 1)].values(), dtype=float), 0
            )

        if np.any(nullity == 0):

            print('No cycles detected, stopping program')
            break

        elif np.any(res != 0.):

            link_robustness, cuttingEdges = calc_link_robustness(
                ref_graph, cut_graph, cyc_nx_base, numeric_res
                )

            robustnessSets.append(link_robustness)
            cutSets.append(cuttingEdges)

        else:

            print('No initial linkages detected, stopping program')
            break

    return np.divide(robustnessSets, nullitySets), cutSets


def calc_link_robustness(ref_graph, cut_graph, cyc_nx_base, numeric_res):
    """
    Sampling randonly through edges and cutting them until the
    intertwined graphs disentangle.

    Args:
        ref_graph (list):\n
            The reference graph, which is setting the frame to be disentangled
            from.
        cut_graph (list):\n
            The graph for which the optimal cut set if searched.
        cyc_nx_base (list):\n
            The initial cycle bases.
        numeric_res (lndarray):\n
            The initial linking number matrix.

    Returns:
        int:\n
            The final cutting number.
        list:\n
            The list of random cuts for the cut_graph.

    """
    link_robustness = -1
    cuttingEdges = []

    # set up dummy graph
    SG = nx.Graph(cut_graph)
    graph_sets = [SG, ref_graph]
    null = itws.calc_nullity(SG)

    while null > 0:

        # cut but do not disconnect
        lk_mat = itws.get_linkage_matrix(numeric_res, cyc_nx_base)
        cycles_to_cut = [j for j, lk in enumerate(lk_mat) if np.any(lk != 0.)]

        # compute linkage of basis vectors
        cyc_id = rd.choice(cycles_to_cut)
        cut_cycle = cyc_nx_base[0][cyc_id]
        eList = list(cut_cycle.edges())
        rd.shuffle(eList)

        for e in eList:

            graph_sets[0].remove_edge(*e)

            if nx.is_connected(graph_sets[0]):

                link_robustness += 1
                cuttingEdges.append(e)

                break

            else:
                graph_sets[0].add_edge(*e)

        cyc_nx_base = [itws.calc_cycle_basis(graph_sets[0]), cyc_nx_base[1]]
        numeric_res = itws.basis_linkage(*cyc_nx_base)

        # check breaks
        null = itws.calc_nullity(graph_sets[0])
        res = np.round(
            np.fromiter(numeric_res[(0, 1)].values(), dtype=float), 0
            )

        if np.all(res == 0):
            break

    return link_robustness, cuttingEdges
