"""Gate 16: bridge the 2-20 um decade and attack discretization directly.

Gate 15 killed the proposed sharp 20-1000 um support threshold under an equal
structural budget.  The remaining loophole is obvious: Gate 13's failed point
coordinate was ~2 um, while Gate 15 started near 20 um.

This gate therefore fixes the support *location* at the driven terminal route and
sweeps physical widths 2-128 um with a much smaller absolute edit budget.  It
then repeats the same experiment after subdividing every morphology edge in two,
so the same physical support is represented by roughly twice as many numerical
compartments.

If the response disappears only because the original morphology file had a
coarse sample interval, original and refined curves should disagree.  If the
response is a physical property of the compiled cable, the curves should agree
when plotted against micrometres even though the compartment counts differ.

The support rule is deliberately no longer forward*return. Gate 15 showed that
selector has not earned credit-assignment status. Here we test only the scale of
operator sensitivity, with a fixed tip-anchored support.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from gate13_real_morphology_length import (
    PointTree,
    SOURCE_NAME,
    choose_distal_tip,
    download_source,
    edge_lengths_um,
    load_dendritic_tree,
    tone_transfer,
    transfer_vector,
)
from gate15_support_response import bounded_score, local_complex_electrotonic_increment


WIDTHS_UM = (2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0)
EDIT_BUDGET_UM = 0.2
OBSERVER_NOISE_FRACTION = 0.01
REFINE_FACTOR = 2


def refine_tree(tree: PointTree, factor: int = REFINE_FACTOR) -> tuple[PointTree, np.ndarray]:
    """Subdivide every physical edge into `factor` serial compartments.

    Returns the refined tree and an old-node -> refined-node map.  Original point
    coordinates remain explicit endpoints, total physical cable length is
    preserved, and radii are linearly interpolated along each old edge.
    """
    if factor < 2:
        raise ValueError("factor must be >= 2")
    pos = [np.asarray(tree.positions_um[0], dtype=float)]
    parents = [-1]
    radii = [float(tree.radii_um[0])]
    stypes = [int(tree.section_types[0])]
    old_to_new = np.zeros(len(tree.parents), dtype=np.int64)
    old_to_new[0] = 0

    for i in range(1, len(tree.parents)):
        p_old = int(tree.parents[i])
        p_new = int(old_to_new[p_old])
        p0 = np.asarray(tree.positions_um[p_old], dtype=float)
        p1 = np.asarray(tree.positions_um[i], dtype=float)
        r0 = float(tree.radii_um[p_old])
        r1 = float(tree.radii_um[i])
        current_parent = p_new
        for k in range(1, factor + 1):
            a = float(k) / float(factor)
            idx = len(pos)
            pos.append((1.0 - a) * p0 + a * p1)
            parents.append(current_parent)
            radii.append((1.0 - a) * r0 + a * r1)
            stypes.append(int(tree.section_types[i]))
            current_parent = idx
        old_to_new[i] = current_parent

    refined = PointTree(
        positions_um=np.asarray(pos, dtype=float),
        parents=np.asarray(parents, dtype=np.int64),
        radii_um=np.asarray(radii, dtype=float),
        section_types=np.asarray(stypes, dtype=np.int64),
    )
    return refined, old_to_new


def tip_anchored_support(path: np.ndarray, lengths_um: np.ndarray, width_um: float) -> np.ndarray:
    """Contiguous terminal support whose cable length is nearest the requested width."""
    path = np.asarray(path, dtype=np.int64)
    nodes = path[1:]  # edge node labels, soma root excluded
    edge = np.asarray(lengths_um, dtype=float)
    if len(nodes) == 0:
        raise ValueError("empty path")
    target = float(width_um)
    cumulative = 0.0
    chosen: list[int] = []
    best: list[int] | None = None
    best_err = float("inf")
    for node in nodes[::-1]:
        # Candidate before adding this edge.
        if chosen:
            err = abs(cumulative - target)
            if err < best_err:
                best = list(chosen)
                best_err = err
        chosen.append(int(node))
        cumulative += float(edge[int(node)])
        err = abs(cumulative - target)
        if err < best_err:
            best = list(chosen)
            best_err = err
        if cumulative >= target and len(chosen) >= 1:
            break
    if best is None:
        best = chosen
    return np.asarray(best[::-1], dtype=np.int64)


def equal_budget_scales(cable_um: float) -> tuple[float, float, float]:
    frac = float(EDIT_BUDGET_UM) / float(cable_um)
    if not 0.0 < frac < 0.5:
        raise ValueError(f"invalid fractional edit {frac}")
    return 1.0 - frac, 1.0, 1.0 + frac


def scaled_support_tones(tree, lengths, tones, tip, support, scale):
    trial = np.asarray(lengths, dtype=float).copy()
    trial[np.asarray(support, dtype=np.int64)] *= float(scale)
    return tone_transfer(tree, trial, np.asarray(tones, dtype=float), int(tip))


def evaluate_width(tree, lengths, path, tip, tones, z0, noise, omega, width_um):
    support = tip_anchored_support(path, lengths, float(width_um))
    cable = float(np.sum(np.asarray(lengths)[support]))
    q0 = bounded_score(z0, noise)
    rows = []
    best = None
    for scale in equal_budget_scales(cable):
        if abs(scale - 1.0) < 1e-15:
            z = np.asarray(z0, dtype=complex)
        else:
            z = scaled_support_tones(tree, lengths, tones, tip, support, scale)
        q = bounded_score(z, noise)
        row = {
            "scale": float(scale),
            "relative_delta_score": float((q - q0) / max(abs(q0), 1e-300)),
            "target_gain_vs_baseline": float(abs(z[2]) / max(abs(z0[2]), 1e-300)),
        }
        rows.append(row)
        if best is None or q > best[0]:
            best = (float(q), row)
    assert best is not None
    gamma_inc = local_complex_electrotonic_increment(tree, lengths, omega)
    return {
        "requested_width_um": float(width_um),
        "support_cable_um": cable,
        "support_nodes": int(len(support)),
        "first_node": int(support[0]),
        "last_node": int(support[-1]),
        "electrotonic_abs_gamma_width": float(np.sum(np.abs(gamma_inc[support]))),
        "gamma_relative": float(best[1]["relative_delta_score"]),
        "best_direction": "shorter" if best[1]["scale"] < 1.0 else ("longer" if best[1]["scale"] > 1.0 else "stay"),
        "best_target_gain_vs_baseline": float(best[1]["target_gain_vs_baseline"]),
        "trials": rows,
    }


def setup_carrier(tree, lengths, tip, omega):
    tones = np.asarray([0.0, 0.35 * omega, omega, 2.5 * omega, 6.0 * omega], dtype=float)
    z0 = tone_transfer(tree, lengths, tones, tip)
    p0 = np.abs(z0) ** 2
    noise = float(OBSERVER_NOISE_FRACTION * p0[2])
    return tones, z0, noise


def best_scalar_residual(a: np.ndarray, b: np.ndarray) -> float:
    scalar = np.vdot(a, b) / max(float(np.vdot(a, a).real), 1e-300)
    return float(np.linalg.norm(b - scalar * a) / max(np.linalg.norm(b), 1e-300))


def run() -> dict:
    out_dir = Path("results")
    source = download_source(out_dir / "_source" / SOURCE_NAME)
    original = load_dendritic_tree(source)
    original_lengths = edge_lengths_um(original)
    tip0, path0, path_len0 = choose_distal_tip(original, original_lengths)

    # Keep the target frequency fixed from the original physical morphology.
    omegas = np.geomspace(0.5, 500.0, 120)
    gains = []
    for w in omegas:
        v, ais = transfer_vector(original, original_lengths, float(w), tip0)
        gains.append(abs(v[ais]))
    peak_i = int(np.argmax(np.asarray(gains)))
    omega_star = float(omegas[peak_i])

    refined, old_to_new = refine_tree(original, REFINE_FACTOR)
    refined_lengths = edge_lengths_um(refined)
    tip1 = int(old_to_new[tip0])
    # Reconstruct the exact refined route to the mapped physical tip.
    rev = []
    j = tip1
    while j != 0:
        rev.append(j)
        j = int(refined.parents[j])
    rev.append(0)
    path1 = np.asarray(rev[::-1], dtype=np.int64)
    path_len1 = float(np.sum(refined_lengths[path1[1:]]))

    tones0, z0, noise0 = setup_carrier(original, original_lengths, tip0, omega_star)
    tones1, z1, noise1 = setup_carrier(refined, refined_lengths, tip1, omega_star)

    baseline = {
        "omega_star_rad_s": omega_star,
        "original_nodes": int(len(original.parents)),
        "refined_nodes": int(len(refined.parents)),
        "original_path_length_um": float(path_len0),
        "refined_path_length_um": float(path_len1),
        "path_length_relative_error": float(abs(path_len1 - path_len0) / max(path_len0, 1e-300)),
        "target_gain_ratio_refined_over_original": float(abs(z1[2]) / max(abs(z0[2]), 1e-300)),
        "five_tone_best_scalar_residual_refined_vs_original": best_scalar_residual(z0, z1),
        "original_bounded_score": float(bounded_score(z0, noise0)),
        "refined_bounded_score": float(bounded_score(z1, noise1)),
    }

    rows0 = [
        evaluate_width(original, original_lengths, path0, tip0, tones0, z0, noise0, omega_star, w)
        for w in WIDTHS_UM
    ]
    rows1 = [
        evaluate_width(refined, refined_lengths, path1, tip1, tones1, z1, noise1, omega_star, w)
        for w in WIDTHS_UM
    ]

    pairs = []
    for a, b in zip(rows0, rows1):
        pairs.append(
            {
                "requested_width_um": float(a["requested_width_um"]),
                "original_actual_width_um": float(a["support_cable_um"]),
                "refined_actual_width_um": float(b["support_cable_um"]),
                "original_support_nodes": int(a["support_nodes"]),
                "refined_support_nodes": int(b["support_nodes"]),
                "original_gamma_relative": float(a["gamma_relative"]),
                "refined_gamma_relative": float(b["gamma_relative"]),
                "gamma_relative_difference": float(b["gamma_relative"] - a["gamma_relative"]),
                "gamma_ratio_refined_over_original": float(
                    b["gamma_relative"] / max(abs(a["gamma_relative"]), 1e-300)
                ),
                "original_electrotonic_width": float(a["electrotonic_abs_gamma_width"]),
                "refined_electrotonic_width": float(b["electrotonic_abs_gamma_width"]),
            }
        )

    gamma0 = np.asarray([r["gamma_relative"] for r in rows0], dtype=float)
    gamma1 = np.asarray([r["gamma_relative"] for r in rows1], dtype=float)
    curve_cos = float(
        np.dot(gamma0, gamma1)
        / max(float(np.linalg.norm(gamma0) * np.linalg.norm(gamma1)), 1e-300)
    )
    rel_curve_error = float(np.linalg.norm(gamma1 - gamma0) / max(np.linalg.norm(gamma0), 1e-300))

    result = {
        "gate": 16,
        "question": (
            "Does equal-budget operator leverage persist from ~2 to 20 um, and does the physical-width response survive a 2x morphology-edge refinement?"
        ),
        "protocol": {
            "widths_um": [float(x) for x in WIDTHS_UM],
            "equal_budget_total_abs_cable_change_um": float(EDIT_BUDGET_UM),
            "observer_noise_fraction": float(OBSERVER_NOISE_FRACTION),
            "support_rule": "tip-anchored contiguous physical cable; no receipt selector",
            "refine_factor": int(REFINE_FACTOR),
        },
        "baseline_discretization_check": baseline,
        "original_curve": rows0,
        "refined_curve": rows1,
        "paired_curve": pairs,
        "curve_agreement": {
            "cosine": curve_cos,
            "relative_l2_error": rel_curve_error,
            "all_original_best_directions": [r["best_direction"] for r in rows0],
            "all_refined_best_directions": [r["best_direction"] for r in rows1],
        },
        "claim_boundary": (
            "Refinement changes the numerical compartmentalization but retains the same piecewise-linear physical morphology and uniform membrane densities. "
            "Agreement would reject one discretization-artifact explanation; it would not establish a biological plasticity scale or prove convergence under arbitrary coarsening/refinement."
        ),
    }
    out_dir.mkdir(exist_ok=True)
    (out_dir / "gate16_microscale_resampling.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    run()
