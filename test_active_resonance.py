import numpy as np

from active_experiment import run
from active_resonance import QuasiActiveDendrite, log_frequency_grid, peak_info


def test_active_gate_relationships_pass():
    r = run()
    assert r["verdict"]["passes"]


def test_passive_and_active_spectra_are_distinct():
    m = QuasiActiveDendrite()
    ws = log_frequency_grid(0.002, 1.0, 300)
    tip = int(m.branch_a[-1])
    a = m.spectrum(ws, tip, m.soma)
    p = m.passive_spectrum(ws, tip, m.soma)
    assert peak_info(ws, p)["omega"] == 0.0
    assert peak_info(ws, a)["omega"] > 0.0
    assert not np.allclose(a, p)


def test_local_quasi_active_edit_is_rank_one_at_fixed_frequency():
    m = QuasiActiveDendrite()
    w = 0.1
    site = int(m.branch_a[len(m.branch_a) // 2])
    dH = m.edited_transfer_matrix(w, site, 0.06) - m.transfer_matrix(w)
    s = np.linalg.svd(dH, compute_uv=False)
    assert s[1] / s[0] < 1e-9


def test_quasi_active_state_is_stable_and_has_complex_poles():
    m = QuasiActiveDendrite()
    poles = m.poles()
    assert np.max(poles.real) < 0.0
    assert np.sum(np.abs(poles.imag) > 1e-8) >= 2
