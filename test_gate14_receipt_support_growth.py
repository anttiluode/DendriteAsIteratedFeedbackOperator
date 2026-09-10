import numpy as np

from gate14_receipt_support_growth import contiguous_support


def test_contiguous_support_keeps_component_containing_peak():
    path = np.asarray([1, 2, 3, 4, 5, 6], dtype=np.int64)
    joint = np.zeros(8, dtype=float)
    joint[path] = [0.1, 0.6, 0.95, 1.0, 0.7, 0.2]
    support = contiguous_support(path, joint, 0.5)
    assert support.tolist() == [2, 3, 4, 5]


def test_tighter_threshold_shrinks_support():
    path = np.asarray([1, 2, 3, 4, 5], dtype=np.int64)
    joint = np.zeros(7, dtype=float)
    joint[path] = [0.2, 0.8, 1.0, 0.85, 0.3]
    loose = contiguous_support(path, joint, 0.5)
    tight = contiguous_support(path, joint, 0.9)
    assert len(tight) < len(loose)
    assert tight.tolist() == [3]
