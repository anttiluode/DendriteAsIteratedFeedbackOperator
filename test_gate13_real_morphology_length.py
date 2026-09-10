import numpy as np

from gate13_real_morphology_length import (
    PointTree,
    build_operator,
    discrete_local_fixed_point,
)


def _tiny_tree():
    return PointTree(
        positions_um=np.asarray([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [20.0, 0.0, 0.0]]),
        parents=np.asarray([-1, 0, 1], dtype=np.int64),
        radii_um=np.asarray([5.0, 1.0, 1.0]),
        section_types=np.asarray([1, 3, 3], dtype=np.int64),
    )


def test_discrete_fixed_point_walks_to_local_maximum():
    values = np.asarray([0.0, 1.0, 3.0, 2.0, 5.0])
    i, history = discrete_local_fixed_point(values, 1)
    assert i == 2
    assert history == [1, 2]


def test_length_edit_changes_membrane_and_axial_operator_together():
    tree = _tiny_tree()
    lengths = np.asarray([0.0, 10.0, 10.0])
    G0, C0, q0, _ = build_operator(tree, lengths)
    G2, C2, q2, _ = build_operator(tree, lengths, edit_node=2, length_scale=2.0)

    # Fixed radius/density: membrane terms at node 2 double.
    assert np.isclose(C2[2] / C0[2], 2.0)
    assert np.isclose(q2[2] / q0[2], 2.0)

    # The direct parent-child off-diagonal is -g_ax; doubling length halves it.
    assert np.isclose(G2[2, 1] / G0[2, 1], 0.5)
