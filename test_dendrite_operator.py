import numpy as np
from dendrite_operator import BranchedCable
from experiment import run


def test_passive_step_is_stable():
    c = BranchedCable()
    rho = np.max(np.abs(np.linalg.eigvals(c.step_operator(0.5))))
    assert rho < 1.0


def test_grounded_fixed_point_matches_iteration():
    c = BranchedCable()
    P = c.step_operator(0.5)
    cue = np.zeros(c.n)
    cue[c.branch_a[-1]] = 1.0
    alpha = 0.08
    analytic = c.grounded_fixed_point(cue, alpha=alpha, dt=0.5)
    x = np.zeros(c.n)
    for _ in range(800):
        x = alpha * cue + (1 - alpha) * (P @ x)
    assert np.linalg.norm(x - analytic) / np.linalg.norm(analytic) < 1e-12


def test_local_edge_resolvent_edit_is_rank_one():
    c = BranchedCable()
    H0 = c.resolvent(0.15)
    i = int(c.branch_b[c.branch_len // 2 - 1])
    j = int(c.branch_b[c.branch_len // 2])
    H1 = c.edited_edge_resolvent(i, j, 0.18, 0.15)
    s = np.linalg.svd(H1 - H0, compute_uv=False)
    assert s[1] / s[0] < 1e-12


def test_receipt_claims_are_reproducible():
    r = run()
    assert r['G1_iterated_mode_selection']['final_abs_cosine_with_slowest_mode'] > 0.99
    assert r['G1_iterated_mode_selection']['final_modal_effective_dimension'] < 1.1
    assert r['G2_grounded_recursion']['analytic_fixed_point_relative_error'] < 1e-12
    assert r['G2_grounded_recursion']['ungrounded_A_vs_B_late_state_cosine'] > 0.99
    assert r['G3_not_a_scalar_weight']['relative_residual_after_best_scalar_fit'] > 0.1
    assert r['G4_local_edit_global_low_rank']['deltaH_effective_rank'] < 1.000001
