import json
from pathlib import Path
import numpy as np

from dendrite_operator import BranchedCable, cosine, effective_rank


def modal_effective_dimension(cable, x):
    lam, Q = cable.modes()
    y = np.sqrt(cable.C) * x
    coeff = Q.T @ y
    power = coeff * coeff
    if power.sum() == 0:
        return 0.0
    p = power / power.sum()
    nz = p > 1e-16
    return float(np.exp(-(p[nz] * np.log(p[nz])).sum()))


def run():
    c = BranchedCable(branch_len=12)
    dt = 0.5
    P = c.step_operator(dt)
    eigP = np.linalg.eigvals(P)

    # G0: passive stability.
    spectral_radius = float(np.max(np.abs(eigP)))

    # G1: Sigh-style autonomous iteration selects the slowest surviving mode.
    x0 = np.zeros(c.n)
    x0[c.branch_a[-1]] = 1.0 / c.C[c.branch_a[-1]]
    lam, Q = c.modes()
    slow_y = Q[:, 0]

    x = x0.copy()
    modal_dim_initial = modal_effective_dimension(c, x)
    traces = []
    for t in range(601):
        if t in (0, 20, 80, 240, 600):
            traces.append({
                "step": t,
                "norm": float(np.linalg.norm(x)),
                "slow_mode_abs_cosine": abs(cosine(np.sqrt(c.C) * x, slow_y)),
                "modal_effective_dimension": modal_effective_dimension(c, x),
            })
        x = P @ x
    modal_dim_final = traces[-1]["modal_effective_dimension"]
    slow_cos_final = traces[-1]["slow_mode_abs_cosine"]

    # G2: grounded recursion has analytic fixed point and stays cue-specific.
    alpha = 0.08
    cue_a = np.zeros(c.n); cue_a[c.branch_a[-1]] = 1.0
    cue_b = np.zeros(c.n); cue_b[c.branch_b[-1]] = 1.0
    fp_a = c.grounded_fixed_point(cue_a, alpha=alpha, dt=dt)
    fp_b = c.grounded_fixed_point(cue_b, alpha=alpha, dt=dt)

    xa = np.zeros(c.n)
    for _ in range(800):
        xa = alpha * cue_a + (1 - alpha) * (P @ xa)
    fixed_point_error = float(np.linalg.norm(xa - fp_a) / np.linalg.norm(fp_a))
    grounded_cue_cosine = cosine(fp_a, fp_b)

    xa = cue_a.copy(); xb = cue_b.copy()
    for _ in range(1000):
        xa = P @ xa; xb = P @ xb
    ungrounded_late_cue_cosine = cosine(np.sqrt(c.C) * xa, np.sqrt(c.C) * xb)

    # G3: equal injected charge on two branches cannot be represented by one scalar kernel.
    ka = c.impulse_kernel(c.branch_a[-1], steps=120, dt=dt)
    kb = c.impulse_kernel(c.branch_b[-1], steps=120, dt=dt)
    beta = float(np.dot(ka, kb) / np.dot(ka, ka))
    scalar_residual = float(np.linalg.norm(kb - beta * ka) / np.linalg.norm(kb))
    kernel_cosine = cosine(ka, kb)
    peak_a = int(np.argmax(np.abs(ka)))
    peak_b = int(np.argmax(np.abs(kb)))

    # G4: a local edge conductance edit creates a dense but rank-one resolvent update.
    H0 = c.resolvent(s=0.15)
    i = int(c.branch_b[c.branch_len // 2 - 1])
    j = int(c.branch_b[c.branch_len // 2])
    H1 = c.edited_edge_resolvent(i, j, delta_g=0.18, s=0.15)
    dH = H1 - H0
    svals = np.linalg.svd(dH, compute_uv=False)
    rank1_ratio = float(svals[0] / svals.sum())
    erank = effective_rank(svals)
    dense_fraction = float(np.mean(np.abs(dH) > 0.02 * np.max(np.abs(dH))))

    receipt = {
        "model": {"compartments": c.n, "branch_len": c.branch_len, "dt": dt},
        "G0_stability": {
            "step_spectral_radius": spectral_radius,
            "passes": spectral_radius < 1.0,
        },
        "G1_iterated_mode_selection": {
            "initial_modal_effective_dimension": modal_dim_initial,
            "final_modal_effective_dimension": modal_dim_final,
            "final_abs_cosine_with_slowest_mode": slow_cos_final,
            "snapshots": traces,
        },
        "G2_grounded_recursion": {
            "alpha": alpha,
            "analytic_fixed_point_relative_error": fixed_point_error,
            "grounded_A_vs_B_fixed_point_cosine": grounded_cue_cosine,
            "ungrounded_A_vs_B_late_state_cosine": ungrounded_late_cue_cosine,
        },
        "G3_not_a_scalar_weight": {
            "branch_kernel_cosine": kernel_cosine,
            "best_scalar_map": beta,
            "relative_residual_after_best_scalar_fit": scalar_residual,
            "peak_step_A": peak_a,
            "peak_step_B": peak_b,
            "peak_value_A": float(ka[peak_a]),
            "peak_value_B": float(kb[peak_b]),
        },
        "G4_local_edit_global_low_rank": {
            "edited_edge": [i, j],
            "delta_conductance": 0.18,
            "deltaH_effective_rank": erank,
            "top_singular_value_fraction": rank1_ratio,
            "dense_fraction_above_2pct_peak": dense_fraction,
            "singular_values_first5": [float(v) for v in svals[:5]],
        },
    }
    return receipt


if __name__ == "__main__":
    r = run()
    print(json.dumps(r, indent=2))
    out = Path("results")
    out.mkdir(exist_ok=True)
    (out / "receipt.json").write_text(json.dumps(r, indent=2) + "\n")
