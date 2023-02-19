# @Author:  Felix Kramer
# @Date:   2021-09-25T14:20:58+02:00
# @Email:  kramer@mpi-cbg.de
# @Project: go-with-the-flow
# @Last modified by:   felix
# @Last modified time: 2022-07-28T18:31:12+02:00
# @License: MIT
import networkx as nx
import numpy as np
import scipy
import intertwined.sampling as itws
import intertwined.edgePriority as itwe


def getEdgeLinkageOperator(graph_sets):
    """
    Compute linking  matrix and edge priority operator for a given set of
    intertwined graphs.

    Args:
        graph_Sets (list): \n
            A list object containing the graphs to be tested.
    Returns:
        ndarray: \n
            The edge priority operator.
        ndarray: \n
            The linkage matrix.

    """

    # cyc_nx_base = [itws.calc_cycle_minimum_basis(G) for G in graph_sets]
    cyc_nx_base = [itws.calc_cycle_basis(G) for G in graph_sets]
    cyc_nx_baseRe = reformat_basis(cyc_nx_base, graph_sets)

    numeric_res = itws.calc_basis_linkage(*cyc_nx_baseRe)
    lk_mat = itws.extract_linkage_matrix(numeric_res)

    graph_matrices = itws.get_basis_matrices(graph_sets, cyc_nx_base)
    cyc_mat_dir = [cm for cm in graph_matrices[-2]]
    cyc_mat_inv = [scipy.linalg.pinv(cm) for cm in cyc_mat_dir]

    P = np.dot(np.dot(cyc_mat_inv[0].T, lk_mat), cyc_mat_inv[1])

    return P, lk_mat


def reformat_basis(cyc_nx_bases, ref_sets):
    """
    ONLY FOR PREVIOSLY COARSE_GRAINED GRAPHS.

    Handling exception and basis reformating in the case of coarse grained
    graphs (see module spatial-network-analysis/span). The goal is to extend
    the formerly toplogically coarse grained cycles and cast them onto their
    full, original embeddings to ensure error-free linkage number calulation
    via the Gauss map.

    Args:
        cyc_nx_bases (list): \n
            A list of cycle bases to be reformated (the coarse grained graph
            objects.)
        ref_sets (list): \n
            A list of cycle bases to be checked in comparison (non-coarse
            grained graph objects with spatial embedding information and full
            nodal structure)
    Returns:
        list \n
            A list of reformated cycle bases.

    """
    cyc_nx_basesRe = [[], []]

    for i, base in enumerate(cyc_nx_bases):

        cm = ref_sets[i].graph['chordsMerged']
        cmPos = ref_sets[i].graph['chordsPos']

        for cycle in base:

            new_cycle = nx.Graph(cycle)
            deg = [n for n in cycle.nodes() if ref_sets[i].degree(n) >= 2]

            allDeg = deg+[deg[0]]
            for j, n in enumerate(allDeg[:-1]):

                if not cycle.has_edge(n, deg[i+1]):

                    subPath = nx.shortest_path(
                        cycle, source=n, target=deg[i+1]
                        )
                    key, expand = getKeyExpand(subPath, cm)

                    if expand:

                        for m in subPath[1:-1]:
                            new_cycle.remove_node(m)

                        nx.add_path(new_cycle, cm[key])
                        nx.set_node_attributes(new_cycle, cmPos[key], 'pos')

            cyc_nx_basesRe[i].append(new_cycle)

    return cyc_nx_basesRe


def getKeyExpand(subPath, cm):
    """
    Test whether a path may be expanded and return path information.

    Args:
        subPath (list): \n
            A sequence of edges to be tested.
        cm (dict): \n
            A dictionary holding expandable path informations
    Returns:
        tuple: \n
            The reformatted subpath object.
        bool: \n
            Boolean, True or False depending on whether subPath may be
            expanded.

    """
    key = []
    expand = False

    if tuple(subPath) in cm.keys():

        key = subPath
        expand = True

    elif tuple(subPath[::-1]) in cm.keys():

        key = subPath[::-1]
        expand = True

    return tuple(key), expand


def cuttingEdgeAlgorithm(cut_graph, ref_graph):
    """
    Compute the optimal cuts set for a graph in an intertwined structure.

    Args:
        cut_graph (list):\n
            The graph for which the optimal cut set if seached.
        ref_graph (list):\n
            The reference graph, which is set the frame to be disentangled.

    Returns:
        list:\n
            The list of optimal cut set for the cut_graph.

    """
    cuts = []
    G1 = nx.Graph(cut_graph)
    nullity = itws.calc_nullity(G1)

    linked = True
    while linked:

        linkage_components = []
        cut_candidates = {}
        cc = [G1.subgraph(c).copy() for c in nx.connected_components(G1)]

        for i, g in enumerate(cc):

            P, lk_mat = getEdgeLinkageOperator([g, ref_graph])

            if not np.all(lk_mat == 0.):
                linkage_components.append(P)

        if len(linkage_components) == 0:
            linked = False
            break

        for P in linkage_components:
            I = np.diagonal(np.dot(P, P.T))
            maxPriority = np.amax(I)
            idx = np.argmax(I)
            cut_candidates[idx] = maxPriority

        idx = max(cut_candidates, key=cut_candidates.get)
        cuts.append(itwe.find_edge_byIndex(G1, idx))

        G1.remove_edge(*cuts[-1])
        if len(cuts) == nullity:
            linked = False
            break

    return cuts
