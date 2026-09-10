import numpy as np

from dendrite_operator import BranchedCable
from active_dendrite import frequency_response, resonant_state_matrix
from gate5_gate6 import run


def test_delayed_restorative_current_creates_resonance():
    c = BranchedCable(branch_len=12)
    A = resonant_state_matrix(
        c, c.branch_b[c.branch_len // 2 :], g_r=0.2, tau=20.0
    )
    omega = np.concatenate([[0.0], np.logspace(-3.0, 0.5, 160)])
    H = frequency_response(c, A, int(c.branch_b[-1]), c.soma, omega)
    mag = np.abs(H)
    assert np.max(np.linalg.eigvals(A).real) < 0.0
    assert mag.max() > 2.0 * mag[0]
    assert omega[np.argmax(mag)] > 0.0


def test_nmda_gate_is_timing_and_state_dependent():
    r = run()
    g6 = r["G6_nmda_state_dependent_operator"]
    assert g6["matched_over_separated"] > 1.5
    assert g6["same_gate_state_vs_zero_deltaJ_rank"] == 2
    assert g6["state_over_zero_probe_transfer"] > 1.05
    assert g6["locally_stable_in_this_run"]
