"""Gate 14: promote the return receipt from one sampled point to a cable support.

Gate 13 transferred the natural-length rule to a real human dendritic morphology
and found a scale problem: the forward x reciprocal-return maximum landed on a
~2 um morphology sampling interval, only ~0.16% of a 1.24 mm input path.  That
coordinate barely changed the global operator and had no finite grounded length
fixed point.

This gate attacks that failure without retuning the membrane model.

The return receipt is treated as a spatial field along the actual soma-to-input
path.  For several predeclared relative thresholds, we take the contiguous
component containing the receipt maximum and coherently scale *all* cable
segments in that support.  The structural coordinate is therefore a supported
cable region rather than a file-discretization interval.

For every support we ask the same Gate-12 question:

    Does shorter/current/longer converge to an interior fixed point of
        target_power / (off_target_power + observer_noise_floor) ?

No support width is selected by looking at the answer.  The thresholds are
fixed before the length sweeps, and negative/boundary results are retained.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from gate13_real_morphology_length import (
    SOURCE_NAME,
    child_counts,
    choose_distal_tip,
    discrete_local_fixed_point,
    download_source,
    edge_lengths_um,
    load_dendritic_tree,
    tone_transfer,
    transfer_vector,
)


SUPPORT_THRESHOLDS = (0.9, 0.5, 0.1)
NOISE_FRACTIONS = (0.001, 0.01, 0.03)


def contiguous_support(path_nodes: np.ndarray, joint: np.ndarray, threshold: float) -> np.ndarray:
    """Contiguous path component around the peak above threshold*peak."""
    path_nodes = np.asarray(path_nodes, dtype=np.int64)
    vals = np.asarray(joint[path_nodes], dtype=float)
    if len(path_nodes) == 0:
        raise ValueError("empty path")
    peak_k = int(np.argmax(vals))
    cutoff = float(threshold) * float(vals[peak_k])
    lo = peak_k
    hi = peak_k
    while lo > 0 and vals[lo - 1] >= cutoff:
        lo -= 1
    while hi + 1 < len(vals) and vals[hi + 1] >= cutoff:
        hi += 1
    return path_nodes[lo : hi + 1].copy()


def scaled_tone_transfer(tree, lengths, tones, tip, support_nodes, scale):
    trial = np.asarray(lengths, dtype=float).copy()
    trial[np.asarray(support_nodes, dtype=np.int64)] *= float(scale)
    return tone_transfer(tree, trial, tones, tip)


def run() -> dict:
    out_dir = Path("results")
    source = download_source(out_dir / "_source" / SOURCE_NAME)
    tree = load_dendritic_tree(source)
    lengths = edge_lengths_um(tree)
    tip, path, path_length = choose_distal_tip(tree, lengths)

    # Same fixed baseline carrier definition as Gate 13.
    omegas = np.geomspace(0.5, 500.0, 120)
    gains = []
    for w in omegas:
        v, ais = transfer_vector(tree, lengths, float(w), tip)
        gains.append(abs(v[ais]))
    gains = np.asarray(gains, dtype=float)
    peak_i = int(np.argmax(gains))
    omega_star = float(omegas[peak_i])
    tones = np.asarray(
        [0.0, 0.35 * omega_star, omega_star, 2.5 * omega_star, 6.0 * omega_star],
        dtype=float,
    )
    z0 = tone_transfer(tree, lengths, tones, tip)
    p0 = np.abs(z0) ** 2
    baseline_purity = float(p0[2] / np.sum(p0))

    # Receipt field: forward influence from the chosen distal tip multiplied by
    # reciprocal small-signal return from the AIS, restricted to this path.
    forward, ais = transfer_vector(tree, lengths, omega_star, tip)
    reverse, _ = transfer_vector(tree, lengths, omega_star, ais)
    joint = np.abs(forward[: len(tree.parents)]) * np.abs(reverse[: len(tree.parents)])
    eligible_path = path[1:]  # each node represents its parent->node cable segment
    peak_node = int(eligible_path[int(np.argmax(joint[eligible_path]))])
    peak_value = float(joint[peak_node])
    joint_norm = joint / max(peak_value, 1e-300)

    # Cable-integrated receipt mass gives a scale diagnostic; it is not used to
    # choose the threshold supports.
    path_mass = float(np.sum(joint[eligible_path] * lengths[eligible_path]))

    scales = np.geomspace(0.2, 5.0, 61)
    start_i = int(np.argmin(np.abs(np.log(scales))))
    support_results = []

    for threshold in SUPPORT_THRESHOLDS:
        support = contiguous_support(eligible_path, joint_norm, threshold)
        support_cable = float(np.sum(lengths[support]))
        support_mass = float(np.sum(joint[support] * lengths[support]))

        grid = np.asarray(
            [
                scaled_tone_transfer(tree, lengths, tones, tip, support, float(scale))
                for scale in scales
            ],
            dtype=complex,
        )
        powers = np.abs(grid) ** 2
        purities = powers[:, 2] / np.maximum(np.sum(powers, axis=1), 1e-300)

        purity_fp_i, purity_history = discrete_local_fixed_point(purities, start_i)
        purity_broad_i = int(np.argmax(purities))

        grounded_rows = []
        reference_i = None
        for frac in NOISE_FRACTIONS:
            noise = float(frac * p0[2])
            score = powers[:, 2] / np.maximum(
                np.sum(powers, axis=1) - powers[:, 2] + noise,
                1e-300,
            )
            fp_i, history = discrete_local_fixed_point(score, start_i)
            broad_i = int(np.argmax(score))
            grounded_rows.append(
                {
                    "observer_noise_fraction_of_baseline_target_power": float(frac),
                    "fixed_point_scale": float(scales[fp_i]),
                    "fixed_point_is_interior": bool(0 < fp_i < len(scales) - 1),
                    "broad_best_scale": float(scales[broad_i]),
                    "broad_best_is_interior": bool(0 < broad_i < len(scales) - 1),
                    "fixed_point_score": float(score[fp_i]),
                    "target_fraction": float(purities[fp_i]),
                    "target_gain_vs_baseline": float(
                        abs(grid[fp_i, 2]) / max(abs(z0[2]), 1e-300)
                    ),
                    "walk_steps": len(history) - 1,
                }
            )
            if abs(frac - 0.01) < 1e-12:
                reference_i = fp_i

        assert reference_i is not None
        zr = grid[int(reference_i)]
        scalar = np.vdot(z0, zr) / max(float(np.vdot(z0, z0).real), 1e-300)
        residual = float(np.linalg.norm(zr - scalar * z0) / max(np.linalg.norm(zr), 1e-300))
        phase = float(np.angle(zr[2] / z0[2]))

        support_results.append(
            {
                "relative_receipt_threshold": float(threshold),
                "support_nodes": int(len(support)),
                "support_first_node": int(support[0]),
                "support_last_node": int(support[-1]),
                "support_cable_um": support_cable,
                "support_fraction_of_input_path": float(support_cable / path_length),
                "support_fraction_of_path_receipt_mass": float(
                    support_mass / max(path_mass, 1e-300)
                ),
                "purity_only": {
                    "fixed_point_scale": float(scales[purity_fp_i]),
                    "fixed_point_is_interior": bool(0 < purity_fp_i < len(scales) - 1),
                    "broad_best_scale": float(scales[purity_broad_i]),
                    "broad_best_is_interior": bool(0 < purity_broad_i < len(scales) - 1),
                    "target_fraction": float(purities[purity_fp_i]),
                    "target_gain_vs_baseline": float(
                        abs(grid[purity_fp_i, 2]) / max(abs(z0[2]), 1e-300)
                    ),
                    "walk_steps": len(purity_history) - 1,
                },
                "grounded": grounded_rows,
                "reference_1pct_operator_change": {
                    "fixed_point_scale": float(scales[int(reference_i)]),
                    "support_total_cable_after_scale_um": float(
                        support_cable * scales[int(reference_i)]
                    ),
                    "target_transfer_phase_shift_radians": phase,
                    "five_tone_best_scalar_residual": residual,
                },
            }
        )

    result = {
        "gate": 14,
        "question": "Does a receipt-sized coherent cable region recover a finite grounded natural-length fixed point on the real morphology?",
        "real_morphology": {
            "dendritic_nodes_including_soma_root": int(len(tree.parents)),
            "total_dendritic_cable_um": float(np.sum(lengths)),
            "selected_input_tip": int(tip),
            "selected_input_path_nodes": int(len(path)),
            "selected_input_path_length_um": float(path_length),
            "branch_nodes": int(np.sum(child_counts(tree.parents) >= 2)),
        },
        "baseline_operator": {
            "omega_star_rad_s": omega_star,
            "peak_is_interior_to_frequency_sweep": bool(0 < peak_i < len(omegas) - 1),
            "target_fraction": baseline_purity,
            "target_gain": float(abs(z0[2])),
        },
        "receipt_field": {
            "peak_node": peak_node,
            "peak_is_terminal_tip": bool(peak_node == tip),
            "support_rule": "contiguous component on the input path with normalized |forward|*|return| >= threshold",
            "thresholds_fixed_before_sweeps": [float(x) for x in SUPPORT_THRESHOLDS],
            "path_integrated_receipt_mass": path_mass,
        },
        "supports": support_results,
        "interpretation_rule": (
            "An interior result counts only if the local neighbour walk and the broad grid optimum are both interior. "
            "Boundary-seeking results remain failures for that structural coordinate."
        ),
        "sigh_bridge": (
            "Gate 13 showed that a point receipt can name a discretization interval too small to be a meaningful growth coordinate. "
            "Gate 14 asks whether the receipt has a spatial support whose coherent metric change can itself approach a structural fixed point."
        ),
        "claim_boundary": (
            "The morphology is real and the support is selected from the baseline transfer field, but the quasi-active membrane and the growth/search rule are phenomenological. "
            "No literal waveguide wavelength, bAP gradient, or biological structural-plasticity mechanism is claimed."
        ),
    }
    out_dir.mkdir(exist_ok=True)
    (out_dir / "gate14_receipt_support_growth.json").write_text(
        json.dumps(result, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    run()
