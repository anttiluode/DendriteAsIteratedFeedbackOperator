"""Gate 17: geometry modes, curvature, and a fair receipt-projection test.

Gates 15-16 killed the attractive claim that dendritic geometry has a sharp
minimum editable support width.  Local operator sensitivity survives down to
micron-scale supports and survives a 2x morphology refinement.

The surviving question is collective:

    Do many locally sensitive length coordinates organize into a small number
    of consequential geometry-change modes?

This gate works on the same pinned human L2/3 morphology as Gates 13-16.  The
soma-to-most-distal-tip route is parameterized by an orthonormal set of smooth
cosine geometry modes in *log length*.  Equal mode amplitude therefore means
an equal length-weighted RMS fractional geometry perturbation, independent of
how many morphology samples happen to represent the route.

We measure:

1. central-difference geometry gradients for 16 spatial modes;
2. an 8x8 finite-difference structural Hessian and its eigenspectrum;
3. projections of forward activity and an event-conditioned return receipt
   onto the same 16 geometry modes;
4. whether those receipt projections predict the *measured* geometry gradient
   better than forward-only / return-only controls and circular-shift nulls.

Important boundaries:
- bra-ket notation is only linear algebra here; no quantum claim;
- the post-event conductance state is phenomenological, not a bAP model;
- overlap does not equal causal credit unless it predicts useful geometry
  directions against controls;
- for maximizing Q, negative Hessian eigenvalues indicate locally concave
  directions.  Positive eigenvalues are not stable maximizing directions.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from gate13_real_morphology_length import (
    SOURCE_NAME,
    TAU_Q_S,
    G_AIS_LEAK_S,
    build_operator,
    choose_distal_tip,
    download_source,
    edge_lengths_um,
    grounded_score,
    load_dendritic_tree,
    tone_transfer,
    transfer_vector,
)


N_GRAD_MODES = 16
N_HESS_MODES = 8
EPS = 0.01  # 1% length-weighted RMS log-length perturbation per normalized mode
OBSERVER_NOISE_FRACTION = 0.01
EVENT_PROXIMAL_UM = 120.0
EVENT_G_DENSITY_S_PER_CM2 = 5.0e-4
EVENT_AIS_SHUNT_S = 5.0e-8
N_SHIFT_NULL = 63


def weighted_cosine_basis(midpoint_fraction: np.ndarray, weights: np.ndarray, n_modes: int) -> np.ndarray:
    """Length-weighted orthonormal cosine basis on an irregularly sampled path."""
    s = np.asarray(midpoint_fraction, dtype=float)
    w = np.asarray(weights, dtype=float)
    if s.ndim != 1 or w.ndim != 1 or len(s) != len(w):
        raise ValueError("s and weights must be same-length vectors")
    if n_modes < 1 or n_modes > len(s):
        raise ValueError("invalid n_modes")
    if np.any(w <= 0):
        raise ValueError("weights must be positive")
    w = w / np.sum(w)

    out: list[np.ndarray] = []
    for k in range(n_modes):
        v = np.cos(np.pi * float(k) * s)
        for q in out:
            v = v - float(np.sum(w * q * v)) * q
        norm = float(np.sqrt(np.sum(w * v * v)))
        if norm < 1e-12:
            raise ValueError(f"degenerate geometry mode {k}")
        out.append(v / norm)
    return np.asarray(out, dtype=float)


def apply_log_length_modes(
    lengths_um: np.ndarray,
    path_segments: np.ndarray,
    basis: np.ndarray,
    coeffs: np.ndarray,
) -> np.ndarray:
    """Apply dimensionless mode coefficients as multiplicative log-length edits."""
    lengths = np.asarray(lengths_um, dtype=float).copy()
    path_segments = np.asarray(path_segments, dtype=np.int64)
    basis = np.asarray(basis, dtype=float)
    coeffs = np.asarray(coeffs, dtype=float)
    if basis.ndim != 2 or basis.shape[1] != len(path_segments):
        raise ValueError("basis/path mismatch")
    if len(coeffs) > basis.shape[0]:
        raise ValueError("too many coefficients")
    log_change = coeffs @ basis[: len(coeffs)]
    lengths[path_segments] *= np.exp(log_change)
    return lengths


def path_coordinates(lengths_um: np.ndarray, path_segments: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    seg = np.asarray(lengths_um, dtype=float)[np.asarray(path_segments, dtype=np.int64)]
    total = float(np.sum(seg))
    if total <= 0:
        raise ValueError("non-positive path length")
    starts = np.concatenate(([0.0], np.cumsum(seg[:-1])))
    mids = starts + 0.5 * seg
    return mids / total, seg / total, mids


def event_return_vector(tree, lengths_um: np.ndarray, omega: float, source_node: int, path_segments: np.ndarray) -> tuple[np.ndarray, int]:
    """Return transfer through a fixed phenomenological post-event conductance state.

    The post-event operator remains reciprocal at this frozen state.  What
    changes is the conductance state relative to the pre-event operator: a
    proximal 120 um dendritic shunt plus an AIS shunt.  This is deliberately a
    simple AHP-like perturbation, not a claim about exact bAP conductances.
    """
    G, cap, gq, ais = build_operator(tree, np.asarray(lengths_um, dtype=float))
    w = float(omega)
    A = G.astype(complex) + sparse.diags(
        1j * w * cap + gq / (1.0 + 1j * w * TAU_Q_S), format="csr"
    )

    event_diag = np.zeros(A.shape[0], dtype=float)
    path_segments = np.asarray(path_segments, dtype=np.int64)
    cumulative = 0.0
    for node in path_segments:
        L = float(lengths_um[int(node)])
        midpoint = cumulative + 0.5 * L
        cumulative += L
        if midpoint > EVENT_PROXIMAL_UM:
            continue
        radius = max(float(tree.radii_um[int(node)]), 0.05)
        area_cm2 = 2.0 * np.pi * radius * max(L, 1e-9) * 1.0e-8
        event_diag[int(node)] += EVENT_G_DENSITY_S_PER_CM2 * area_cm2
    event_diag[ais] += EVENT_AIS_SHUNT_S
    A = A + sparse.diags(event_diag.astype(complex), format="csr")

    b = np.zeros(A.shape[0], dtype=complex)
    b[int(source_node)] = 1.0
    return np.asarray(spsolve(A.tocsc(), b), dtype=complex), ais


def safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ac = a - np.mean(a)
    bc = b - np.mean(b)
    den = float(np.linalg.norm(ac) * np.linalg.norm(bc))
    if den <= 1e-30:
        return 0.0
    return float(np.dot(ac, bc) / den)


def projection_scores(basis: np.ndarray, weights: np.ndarray, density: np.ndarray) -> np.ndarray:
    basis = np.asarray(basis, dtype=float)
    weights = np.asarray(weights, dtype=float)
    density = np.asarray(density, dtype=float)
    return np.asarray([np.sum(weights * q * density) for q in basis], dtype=float)


def circular_shift_null(
    basis: np.ndarray,
    weights: np.ndarray,
    density: np.ndarray,
    target_gradient: np.ndarray,
    n_shift: int = N_SHIFT_NULL,
) -> dict:
    n = len(density)
    shifts = np.unique(np.linspace(1, max(1, n - 1), int(n_shift), dtype=int))
    vals = []
    for shift in shifts:
        score = projection_scores(basis, weights, np.roll(density, int(shift)))
        vals.append(abs(safe_corr(score, target_gradient)))
    vals = np.asarray(vals, dtype=float)
    actual = abs(safe_corr(projection_scores(basis, weights, density), target_gradient))
    p = float((1.0 + np.sum(vals >= actual)) / (1.0 + len(vals)))
    return {
        "absolute_actual_correlation": actual,
        "circular_shift_empirical_p": p,
        "null_median_abs_correlation": float(np.median(vals)) if len(vals) else 0.0,
        "null_max_abs_correlation": float(np.max(vals)) if len(vals) else 0.0,
        "n_shifts": int(len(vals)),
    }


def effective_rank_from_spectrum(values: np.ndarray) -> float:
    x = np.abs(np.asarray(values, dtype=float))
    s = float(np.sum(x))
    if s <= 1e-30:
        return 0.0
    p = x / s
    p = p[p > 0]
    return float(np.exp(-np.sum(p * np.log(p))))


def energy_rank(values: np.ndarray, fraction: float = 0.95) -> int:
    e = np.sort(np.asarray(values, dtype=float) ** 2)[::-1]
    total = float(np.sum(e))
    if total <= 1e-30:
        return 0
    return int(np.searchsorted(np.cumsum(e) / total, fraction) + 1)


def run() -> dict:
    out_dir = Path("results")
    source = download_source(out_dir / "_source" / SOURCE_NAME)
    tree = load_dendritic_tree(source)
    lengths = edge_lengths_um(tree)
    tip, path, path_length = choose_distal_tip(tree, lengths)
    path_segments = np.asarray(path[1:], dtype=np.int64)

    # Keep the same carrier definition used by Gates 13-16.
    omegas = np.geomspace(0.5, 500.0, 120)
    gains = []
    for omega in omegas:
        v, ais = transfer_vector(tree, lengths, float(omega), tip)
        gains.append(abs(v[ais]))
    gains = np.asarray(gains, dtype=float)
    peak_i = int(np.argmax(gains))
    omega_star = float(omegas[peak_i])
    tones = np.asarray([0.0, 0.35 * omega_star, omega_star, 2.5 * omega_star, 6.0 * omega_star])

    z0 = tone_transfer(tree, lengths, tones, tip)
    p0 = np.abs(z0) ** 2
    noise_power = float(OBSERVER_NOISE_FRACTION * p0[2])
    q0 = float(grounded_score(z0, noise_power))

    s, weights, mids_um = path_coordinates(lengths, path_segments)
    basis = weighted_cosine_basis(s, weights, N_GRAD_MODES)
    gram = (basis * weights[None, :]) @ basis.T
    basis_error = float(np.linalg.norm(gram - np.eye(N_GRAD_MODES)))

    def objective(coeffs: np.ndarray) -> float:
        trial = apply_log_length_modes(lengths, path_segments, basis, np.asarray(coeffs, dtype=float))
        z = tone_transfer(tree, trial, tones, tip)
        return float(grounded_score(z, noise_power))

    plus = np.zeros(N_GRAD_MODES, dtype=float)
    minus = np.zeros(N_GRAD_MODES, dtype=float)
    gradient = np.zeros(N_GRAD_MODES, dtype=float)
    diag_curvature = np.zeros(N_GRAD_MODES, dtype=float)
    for k in range(N_GRAD_MODES):
        c = np.zeros(N_GRAD_MODES, dtype=float)
        c[k] = EPS
        plus[k] = objective(c)
        c[k] = -EPS
        minus[k] = objective(c)
        gradient[k] = (plus[k] - minus[k]) / (2.0 * EPS)
        diag_curvature[k] = (plus[k] - 2.0 * q0 + minus[k]) / (EPS * EPS)

    # Full Hessian only in the first eight smooth modes.  Diagonal terms reuse
    # the already measured central differences above.
    m = N_HESS_MODES
    K = np.zeros((m, m), dtype=float)
    for i in range(m):
        K[i, i] = diag_curvature[i]
    for i in range(m):
        for j in range(i + 1, m):
            vals = {}
            for si in (-1.0, 1.0):
                for sj in (-1.0, 1.0):
                    c = np.zeros(N_GRAD_MODES, dtype=float)
                    c[i] = si * EPS
                    c[j] = sj * EPS
                    vals[(si, sj)] = objective(c)
            kij = (
                vals[(1.0, 1.0)]
                - vals[(1.0, -1.0)]
                - vals[(-1.0, 1.0)]
                + vals[(-1.0, -1.0)]
            ) / (4.0 * EPS * EPS)
            K[i, j] = kij
            K[j, i] = kij

    evals, evecs = np.linalg.eigh(K)
    order = np.argsort(np.abs(evals))[::-1]
    evals_ordered = evals[order]
    evecs_ordered = evecs[:, order]
    projected_gradient = evecs_ordered.T @ gradient[:m]
    stationary_candidates = []
    for rank, (lam, gg) in enumerate(zip(evals_ordered, projected_gradient), start=1):
        a_star = None if abs(lam) < 1e-20 else float(-gg / lam)
        stationary_candidates.append(
            {
                "rank_by_abs_curvature": rank,
                "eigenvalue": float(lam),
                "directional_gradient": float(gg),
                "quadratic_stationary_amplitude": a_star,
                "locally_concave_for_maximization": bool(lam < 0),
                "stationary_point_within_25pct_log_rms": bool(
                    a_star is not None and lam < 0 and abs(a_star) <= 0.25
                ),
                "basis_coefficients": [float(x) for x in evecs_ordered[:, rank - 1]],
            }
        )

    # Activity / return fields at the same baseline target carrier.
    forward, ais = transfer_vector(tree, lengths, omega_star, tip)
    reverse_pre, _ = transfer_vector(tree, lengths, omega_star, ais)
    reverse_post, _ = event_return_vector(tree, lengths, omega_star, ais, path_segments)
    delta_return = reverse_post - reverse_pre

    e = forward[path_segments]
    r0 = reverse_pre[path_segments]
    dr = delta_return[path_segments]
    inner = np.sum(weights * np.conj(dr) * e)
    norm_e = float(np.sum(weights * np.abs(e) ** 2))
    norm_dr = float(np.sum(weights * np.abs(dr) ** 2))
    fidelity = float(abs(inner) ** 2 / max(norm_e * norm_dr, 1e-300))
    residual_vs_reciprocal = float(
        np.linalg.norm(dr - (np.vdot(r0, dr) / max(float(np.vdot(r0, r0).real), 1e-300)) * r0)
        / max(np.linalg.norm(dr), 1e-300)
    )

    densities = {
        "event_phase_overlap": np.real(np.conj(dr) * e),
        "event_magnitude_product": np.abs(dr) * np.abs(e),
        "forward_power": np.abs(e) ** 2,
        "return_residual_power": np.abs(dr) ** 2,
        "reciprocal_phase_overlap": np.real(np.conj(r0) * e),
    }
    predictors = {}
    g_rel = gradient / max(abs(q0), 1e-300)
    for name, density in densities.items():
        scores = projection_scores(basis, weights, density)
        corr = safe_corr(scores, g_rel)
        top = np.argsort(np.abs(g_rel))[::-1][: max(4, N_GRAD_MODES // 4)]
        sign_hits = int(np.sum(np.sign(scores[top]) == np.sign(g_rel[top])))
        predictors[name] = {
            "projection_scores": [float(x) for x in scores],
            "pearson_with_measured_gradient": float(corr),
            "abs_correlation": float(abs(corr)),
            "top_gradient_sign_hits": sign_hits,
            "top_gradient_sign_total": int(len(top)),
            "circular_shift_null": circular_shift_null(basis, weights, density, g_rel),
        }

    grad_energy = gradient * gradient
    low4_fraction = float(np.sum(grad_energy[:4]) / max(float(np.sum(grad_energy)), 1e-300))
    high4_fraction = float(np.sum(grad_energy[-4:]) / max(float(np.sum(grad_energy)), 1e-300))

    result = {
        "gate": 17,
        "question": "Do locally sensitive dendritic length coordinates organize into collective geometry modes, and does an event-conditioned return predict their useful direction?",
        "protocol": {
            "real_morphology_nodes": int(len(tree.parents)),
            "selected_tip": int(tip),
            "selected_path_length_um": float(path_length),
            "omega_star_rad_s": omega_star,
            "observer_noise_fraction_of_baseline_target_power": OBSERVER_NOISE_FRACTION,
            "geometry_parameterization": "length-weighted orthonormal cosine modes in log segment length",
            "gradient_modes": N_GRAD_MODES,
            "hessian_modes": N_HESS_MODES,
            "finite_difference_amplitude": EPS,
            "basis_orthonormality_frobenius_error": basis_error,
            "post_event_state": {
                "proximal_path_extent_um": EVENT_PROXIMAL_UM,
                "extra_conductance_density_S_per_cm2": EVENT_G_DENSITY_S_PER_CM2,
                "AIS_extra_shunt_S": EVENT_AIS_SHUNT_S,
            },
        },
        "baseline": {
            "grounded_score": q0,
            "target_fraction": float(p0[2] / np.sum(p0)),
            "target_gain": float(abs(z0[2])),
        },
        "geometry_gradient": {
            "relative_gradient_per_unit_mode_amplitude": [float(x) for x in g_rel],
            "relative_diagonal_curvature_per_unit2": [float(x / max(abs(q0), 1e-300)) for x in diag_curvature],
            "gradient_energy_fraction_first_4_smooth_modes": low4_fraction,
            "gradient_energy_fraction_last_4_fine_modes": high4_fraction,
            "strongest_mode_by_abs_gradient": int(np.argmax(np.abs(gradient))),
        },
        "structural_hessian_first_8_modes": {
            "relative_matrix": [[float(x / max(abs(q0), 1e-300)) for x in row] for row in K],
            "relative_eigenvalues_by_abs_curvature": [float(x / max(abs(q0), 1e-300)) for x in evals_ordered],
            "negative_eigenvalues": int(np.sum(evals < 0)),
            "positive_eigenvalues": int(np.sum(evals > 0)),
            "effective_rank_abs_spectrum": effective_rank_from_spectrum(evals),
            "frobenius_energy_rank_95pct": energy_rank(evals, 0.95),
            "stationary_candidates": stationary_candidates,
        },
        "receipt_bra_ket_test": {
            "weighted_inner_product_real": float(np.real(inner)),
            "weighted_inner_product_imag": float(np.imag(inner)),
            "normalized_fidelity": fidelity,
            "event_return_residual_non_scalar_vs_reciprocal_return": residual_vs_reciprocal,
            "predictors": predictors,
        },
        "interpretation_rule": (
            "A low-dimensional geometry claim requires concentrated gradient/curvature spectra rather than merely nonzero local sensitivity. "
            "A receipt-credit claim requires event-conditioned projections to predict measured geometry gradients better than forward-only/return-only controls and spatial-shift nulls."
        ),
        "claim_boundary": (
            "The morphology is real.  Membrane kinetics, post-event conductance state, observer objective, and geometry modes are phenomenological. "
            "Bra-ket notation denotes ordinary complex vector spaces; no quantum mechanism is claimed."
        ),
    }
    out_dir.mkdir(exist_ok=True)
    (out_dir / "gate17_geometry_modes.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    run()
