# @Author: Felix Kramer <kramer>
# @Date:   26-07-2022
# @Email:  felixuwekramer@proton.me
# @Last modified by:   felix
# @Last modified time: 2022-07-28T10:52:38+02:00
import numpy as np
import snarlpy.sampling as sps
import auxFunc as af


# compute linking numbers for unlink, simple link sand double link and assert
# correctness for toy systems
def test_linkingNumbers():

    results = [
        [[0, 0]],
        [[-1, 0]],
        [[-1, -1]],
        [[1, -1]],
    ]

    for i, constr in enumerate(af.constr_gen):

        K = constr()
        graphSets = [k.G for k in K.layer]
        cyc_nx_base, lk_mat = sps.calc_basisIntertwinedness(graphSets)
        assert np.array_equal(results[i], lk_mat)
