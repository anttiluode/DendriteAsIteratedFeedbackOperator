import numpy as np

from gate15_support_response import (
    candidate_windows,
    equal_budget_scales,
    overlap_fraction,
    select_supports,
)


def test_equal_budget_scales_fix_total_change():
    shorter, current, longer = equal_budget_scales(320.0, budget_um=4.0)
    assert current == 1.0
    assert np.isclose((1.0 - shorter) * 320.0, 4.0)
    assert np.isclose((longer - 1.0) * 320.0, 4.0)


def test_candidate_windows_track_requested_width():
    path = np.arange(1, 101, dtype=np.int64)
    lengths = np.zeros(101, dtype=float)
    lengths[1:] = 2.0
    rows = candidate_windows(path, lengths, 40.0)
    assert rows
    assert min(abs(w - 40.0) for _, _, w in rows) < 1e-12


def test_selectors_return_contiguous_equal_scale_supports():
    path = np.arange(1, 101, dtype=np.int64)
    lengths = np.zeros(101, dtype=float)
    lengths[1:] = 2.0
    f = np.zeros(101, dtype=float)
    r = np.zeros(101, dtype=float)
    f[70:101] = np.linspace(0.1, 1.0, 31)
    r[50:101] = np.linspace(0.2, 1.0, 51)
    supports = select_supports(path, lengths, f, r, 40.0, np.random.default_rng(3))
    assert set(supports) == {"forward", "reverse", "joint", "random", "distance_matched"}
    widths = []
    for s in supports.values():
        idx = np.searchsorted(path, s)
        assert np.all(np.diff(idx) == 1)
        widths.append(np.sum(lengths[s]))
    assert max(widths) - min(widths) <= 4.0


def test_overlap_fraction_is_cable_weighted():
    lengths = np.zeros(10, dtype=float)
    lengths[1:] = 2.0
    a = np.array([1, 2, 3, 4])
    b = np.array([3, 4, 5, 6])
    assert np.isclose(overlap_fraction(a, b, lengths), 0.5)
