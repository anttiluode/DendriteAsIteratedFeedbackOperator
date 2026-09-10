import numpy as np

from gate13_real_morphology_length import PointTree, edge_lengths_um
from gate16_microscale_resampling import refine_tree, tip_anchored_support


def tiny_tree():
    return PointTree(
        positions_um=np.array([[0,0,0],[2,0,0],[5,0,0],[9,0,0]], dtype=float),
        parents=np.array([-1,0,1,2], dtype=np.int64),
        radii_um=np.array([2.0,1.0,0.8,0.6], dtype=float),
        section_types=np.array([1,3,3,3], dtype=np.int64),
    )


def test_refinement_preserves_endpoints_and_total_cable():
    tree = tiny_tree()
    refined, old_to_new = refine_tree(tree, 2)
    l0 = edge_lengths_um(tree)
    l1 = edge_lengths_um(refined)
    assert len(refined.parents) == 1 + 2 * (len(tree.parents) - 1)
    assert np.isclose(np.sum(l0), np.sum(l1))
    for i in range(len(tree.parents)):
        assert np.allclose(refined.positions_um[old_to_new[i]], tree.positions_um[i])


def test_tip_support_tracks_physical_width_after_refinement():
    tree = tiny_tree()
    lengths = edge_lengths_um(tree)
    path = np.array([0,1,2,3], dtype=np.int64)
    s0 = tip_anchored_support(path, lengths, 4.0)
    refined, old_to_new = refine_tree(tree, 2)
    lr = edge_lengths_um(refined)
    tip = int(old_to_new[3])
    rev=[]; j=tip
    while j != 0:
        rev.append(j); j=int(refined.parents[j])
    rev.append(0)
    pr=np.array(rev[::-1], dtype=np.int64)
    s1 = tip_anchored_support(pr, lr, 4.0)
    assert abs(np.sum(lengths[s0]) - 4.0) <= 2.0
    assert abs(np.sum(lr[s1]) - 4.0) <= 1.0
    assert len(s1) >= len(s0)
