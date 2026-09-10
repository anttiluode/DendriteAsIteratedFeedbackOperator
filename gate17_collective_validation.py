"""Gate 17 validation: does the measured low-mode quadratic point survive replay?

The main Gate 17 experiment measures the first eight-mode geometry Hessian around
baseline.  This sidecar does not refit that Hessian.  It reads the frozen-in-run
Gate 17 measurement, forms the Newton candidate -K^{-1} g in the same eight-mode
subspace, then directly re-evaluates the nonlinear cable model along that
collective direction and probes the local gradient at the best counted point.

This distinguishes a suggestive local curvature spectrum from an actually useful
collective structural fixed point.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from gate13_real_morphology_length import (
    SOURCE_NAME,
    choose_distal_tip,
    download_source,
    edge_lengths_um,
    grounded_score,
    load_dendritic_tree,
    tone_transfer,
)
from gate17_geometry_modes import (
    EPS,
    N_GRAD_MODES,
    N_HESS_MODES,
    OBSERVER_NOISE_FRACTION,
    apply_log_length_modes,
    path_coordinates,
    weighted_cosine_basis,
)


LINE_FRACTIONS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5)


def run() -> dict:
    out_dir = Path("results")
    measured_path = out_dir / "gate17_geometry_modes.json"
    if not measured_path.exists():
        raise FileNotFoundError("run gate17_geometry_modes.py first")
    measured = json.loads(measured_path.read_text(encoding="utf-8"))

    source = download_source(out_dir / "_source" / SOURCE_NAME)
    tree = load_dendritic_tree(source)
    lengths = edge_lengths_um(tree)
    tip, path, path_length = choose_distal_tip(tree, lengths)
    path_segments = np.asarray(path[1:], dtype=np.int64)

    omega_star = float(measured["protocol"]["omega_star_rad_s"])
    tones = np.asarray([0.0, 0.35 * omega_star, omega_star, 2.5 * omega_star, 6.0 * omega_star])
    z0 = tone_transfer(tree, lengths, tones, tip)
    p0 = np.abs(z0) ** 2
    noise_power = float(OBSERVER_NOISE_FRACTION * p0[2])
    q0 = float(grounded_score(z0, noise_power))

    s, weights, _ = path_coordinates(lengths, path_segments)
    basis = weighted_cosine_basis(s, weights, N_GRAD_MODES)

    K_rel = np.asarray(measured["structural_hessian_first_8_modes"]["relative_matrix"], dtype=float)
    g_rel = np.asarray(measured["geometry_gradient"]["relative_gradient_per_unit_mode_amplitude"], dtype=float)[:N_HESS_MODES]
    if K_rel.shape != (N_HESS_MODES, N_HESS_MODES):
        raise ValueError("unexpected Hessian shape")

    # Q0 cancels because both K and g are stored relative to Q0.
    newton = -np.linalg.solve(K_rel, g_rel)

    def objective(coeffs8: np.ndarray) -> float:
        coeffs = np.zeros(N_GRAD_MODES, dtype=float)
        coeffs[:N_HESS_MODES] = np.asarray(coeffs8, dtype=float)
        trial = apply_log_length_modes(lengths, path_segments, basis, coeffs)
        z = tone_transfer(tree, trial, tones, tip)
        return float(grounded_score(z, noise_power))

    line = []
    for frac in LINE_FRACTIONS:
        q = objective(float(frac) * newton)
        line.append(
            {
                "fraction_of_newton_candidate": float(frac),
                "relative_score_change": float(q / q0 - 1.0),
                "score": q,
            }
        )
    best = max(line, key=lambda row: row["score"])
    best_frac = float(best["fraction_of_newton_candidate"])
    candidate = best_frac * newton
    q_candidate = float(best["score"])

    # Re-measure first derivatives at the candidate.  If the candidate is a
    # genuine local fixed point in this subspace, these should shrink relative
    # to the baseline gradient rather than merely improving one line direction.
    candidate_gradient = np.zeros(N_HESS_MODES, dtype=float)
    for k in range(N_HESS_MODES):
        plus = candidate.copy()
        minus = candidate.copy()
        plus[k] += EPS
        minus[k] -= EPS
        candidate_gradient[k] = (objective(plus) - objective(minus)) / (2.0 * EPS * q0)

    baseline_grad_norm = float(np.linalg.norm(g_rel))
    candidate_grad_norm = float(np.linalg.norm(candidate_gradient))
    result = {
        "gate": "17-validation",
        "question": "Does the eight-mode quadratic geometry candidate survive direct nonlinear replay and approach a collective local stationary point?",
        "selected_path_length_um": float(path_length),
        "newton_candidate_coefficients": [float(x) for x in newton],
        "newton_candidate_l2_norm": float(np.linalg.norm(newton)),
        "newton_candidate_max_abs_coefficient": float(np.max(np.abs(newton))),
        "counted_line_scan": line,
        "best_line_fraction": best_frac,
        "best_candidate_coefficients": [float(x) for x in candidate],
        "best_relative_score_change": float(q_candidate / q0 - 1.0),
        "baseline_relative_gradient_l2": baseline_grad_norm,
        "candidate_relative_gradient_components": [float(x) for x in candidate_gradient],
        "candidate_relative_gradient_l2": candidate_grad_norm,
        "gradient_norm_ratio_candidate_over_baseline": float(
            candidate_grad_norm / max(baseline_grad_norm, 1e-300)
        ),
        "candidate_is_more_stationary_than_baseline": bool(candidate_grad_norm < baseline_grad_norm),
        "claim_boundary": (
            "This validates only an eight-mode smooth log-length subspace around one phenomenological objective on one real morphology. "
            "It is not a biological growth law or proof that the full morphology is at, or moves toward, a natural shape."
        ),
    }
    (out_dir / "gate17_collective_validation.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    run()
