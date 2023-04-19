# @Author: Felix Kramer <kramer>
# @Date:   26-07-2022
# @Email:  felixuwekramer@proton.me
# @Last modified by:   felix
# @Last modified time: 2022-07-30T09:59:48+02:00
import numpy as np
import snarlpy.edgePriority as spe
import auxFunc as af
# compute linking numbers for unlink, simple link sand double link and assert
# correctness for toy systems


def test_edgeCut():

    results = [
        [],
        [(0, 1)],
        [(0, 1), (1, 2)],
        [(1, 4)],
    ]

    for i, constr in enumerate(af.constr_gen):

        K = constr()
        print(constr)
        graphSets = [k.G for k in K.layer]
        cut_list = spe.cuttingEdgeAlgorithm(*graphSets[::-1])

        print(cut_list)
        assert np.array_equal(results[i], cut_list)
