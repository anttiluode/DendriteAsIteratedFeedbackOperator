import numpy as np

from gate17_geometry_modes import (
    apply_log_length_modes,
    energy_rank,
    projection_scores,
    safe_corr,
    weighted_cosine_basis,
)


def test_weighted_cosine_basis_is_orthonormal_on_irregular_path():
    seg = np.asarray([1.0, 2.0, 4.0, 3.0, 5.0, 2.5, 1.5, 6.0])
    weights = seg / seg.sum()
    starts = np.concatenate(([0.0], np.cumsum(seg[:-1])))
    mids = (starts + 0.5 * seg) / seg.sum()
    basis = weighted_cosine_basis(mids, weights, 6)
    gram = (basis * weights[None, :]) @ basis.T
    assert np.allclose(gram, np.eye(6), atol=1e-12)


def test_equal_mode_amplitude_has_equal_weighted_rms_log_change():
    lengths = np.asarray([0.0, 1.0, 2.0, 4.0, 3.0, 5.0, 2.5, 1.5, 6.0])
    path = np.arange(1, len(lengths), dtype=np.int64)
    seg = lengths[path]
    weights = seg / seg.sum()
    starts = np.concatenate(([0.0], np.cumsum(seg[:-1])))
    mids = (starts + 0.5 * seg) / seg.sum()
    basis = weighted_cosine_basis(mids, weights, 5)
    eps = 0.013
    for k in range(5):
        coeffs = np.zeros(5)
        coeffs[k] = eps
        changed = apply_log_length_modes(lengths, path, basis, coeffs)
        log_change = np.log(changed[path] / lengths[path])
        rms = np.sqrt(np.sum(weights * log_change * log_change))
        assert np.isclose(rms, eps, atol=1e-12)
        assert np.all(changed[path] > 0)


def test_projection_and_correlation_helpers():
    weights = np.asarray([0.25, 0.25, 0.25, 0.25])
    basis = np.asarray([[1, 1, 1, 1], [1, 1, -1, -1]], dtype=float)
    density = np.asarray([2.0, 2.0, -1.0, -1.0])
    scores = projection_scores(basis, weights, density)
    assert np.allclose(scores, [0.5, 1.5])
    assert safe_corr(scores, np.asarray([1.0, 3.0])) > 0.999999


def test_energy_rank_reports_concentrated_spectrum():
    assert energy_rank(np.asarray([10.0, 1.0, 0.1, 0.01]), 0.95) == 1
