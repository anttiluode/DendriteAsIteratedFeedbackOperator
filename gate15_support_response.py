"""Gate 15: measure the spatial support scale at which geometry gains operator leverage.

Gate 13 showed that one ~2 um morphology sample interval was too microscopic.
Gate 14 showed that coherently editing hundreds of micrometres of the same real
soma-to-tip route can produce interior bounded-observer optima.  Gate 15 attacks
the obvious confound: larger supports also edit more cable.

This gate therefore fixes an *absolute structural budget* B.  For a support of
physical cable width W, the two non-zero trial edits are chosen so that

    sum_j |Delta L_j| = B

for every selector and every W.  A larger support receives a proportionally
smaller fractional metric change.  The measured support-response curve is

    Gamma_B(W) = max_{+/- equal-budget edit} [Q_sigma(after) - Q_sigma(before)]

with the unedited state included, so a support that cannot improve the bounded
readout has Gamma_B(W)=0 rather than being forced to move.

Five equal-width support selectors are compared on the pinned human L2/3 cell:

    forward          integrated |H_{x,tip}|
    reverse          integrated |H_{x,AIS}|
    joint            integrated |forward| |reverse|
    random           deterministic random contiguous support
    distance_matched same gross path-position as the joint support, but chosen
                     to have low joint score and little overlap when possible

The same joint-selector width sweep is repeated at 0.5, 1 and 2 times the
baseline preferred frequency.  Every support is also reported in a local
frequency-dependent electrotonic coordinate derived from the same membrane
parameters.  That coordinate is a diagnostic, not a claim that a branched
active dendrite is an exact uniform transmission line.

The gate is deliberately allowed to fail.  If Gamma_B(W) is flat, smooth, or
largest at the smallest support, there is no earned "coherence threshold".
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from gate13_real_morphology_length import (
    CM_F_PER_CM2,
    GLEAK_S_PER_CM2,
    GQ_S_PER_CM2,
    RA_OHM_CM,
    SOURCE_NAME,
    TAU_Q_S,
    choose_distal_tip,
    download_source,
    edge_lengths_um,
    load_dendritic_tree,
    tone_transfer,
    transfer_vector,
)


WIDTHS_UM = (20.0, 40.0, 80.0, 160.0, 320.0, 640.0, 1000.0)
EDIT_BUDGET_UM = 4.0
OBSERVER_NOISE_FRACTION = 0.01
FREQUENCY_FACTORS = (0.5, 1.0, 2.0)
RANDOM_SEED = 15015
SELECTORS = ("forward", "reverse", "joint", "random", "distance_matched")


def path_edge_midpoints_um(path_nodes: np.ndarray, lengths_um: np.ndarray) -> np.ndarray:
    """Distance from soma to the midpoint of every path edge."""
    nodes = np.asarray(path_nodes, dtype=np.int64)
    lens = np.asarray(lengths_um, dtype=float)[nodes]
    starts = np.concatenate(([0.0], np.cumsum(lens[:-1])))
    return starts + 0.5 * lens


def candidate_windows(path_nodes: np.ndarray, lengths_um: np.ndarray, target_width_um: float):
    """Near-equal-width contiguous supports spanning a target physical cable length.

    For every start edge, retain the window ending immediately before and after
    the cumulative width crosses the target.  Then keep windows within one local
    sample interval (or 2%, whichever is larger) of the best width mismatch.
    """
    nodes = np.asarray(path_nodes, dtype=np.int64)
    edge = np.asarray(lengths_um, dtype=float)[nodes]
    if len(nodes) == 0:
        raise ValueError("empty path")
    if target_width_um <= 0:
        raise ValueError("target width must be positive")
    prefix = np.concatenate(([0.0], np.cumsum(edge)))
    raw: list[tuple[int, int, float]] = []
    for i in range(len(nodes)):
        goal = prefix[i] + float(target_width_um)
        j = int(np.searchsorted(prefix, goal, side="left"))
        for jj in (j - 1, j):
            jj = max(i + 1, min(jj, len(nodes)))
            width = float(prefix[jj] - prefix[i])
            raw.append((i, jj, width))
    # Deduplicate identical start/end pairs.
    uniq: dict[tuple[int, int], float] = {}
    for i, j, w in raw:
        uniq[(i, j)] = w
    rows = [(i, j, w) for (i, j), w in uniq.items()]
    mismatch = np.asarray([abs(w - target_width_um) for _, _, w in rows], dtype=float)
    best = float(np.min(mismatch))
    tol = max(0.02 * float(target_width_um), float(np.median(edge)))
    kept = [row for row, d in zip(rows, mismatch) if d <= best + tol]
    return kept


def integrated_window_score(
    path_nodes: np.ndarray,
    lengths_um: np.ndarray,
    field: np.ndarray,
    window: tuple[int, int, float],
) -> float:
    i, j, _ = window
    nodes = np.asarray(path_nodes, dtype=np.int64)[i:j]
    return float(np.sum(np.asarray(field, dtype=float)[nodes] * np.asarray(lengths_um)[nodes]))


def support_from_window(path_nodes: np.ndarray, window: tuple[int, int, float]) -> np.ndarray:
    i, j, _ = window
    return np.asarray(path_nodes, dtype=np.int64)[i:j].copy()


def overlap_fraction(a: np.ndarray, b: np.ndarray, lengths_um: np.ndarray) -> float:
    aset = set(int(x) for x in np.asarray(a, dtype=np.int64))
    bset = set(int(x) for x in np.asarray(b, dtype=np.int64))
    inter = np.asarray(sorted(aset & bset), dtype=np.int64)
    if len(inter) == 0:
        return 0.0
    denom = min(
        float(np.sum(np.asarray(lengths_um)[np.asarray(list(aset), dtype=np.int64)])),
        float(np.sum(np.asarray(lengths_um)[np.asarray(list(bset), dtype=np.int64)])),
    )
    return float(np.sum(np.asarray(lengths_um)[inter]) / max(denom, 1e-300))


def select_supports(
    path_nodes: np.ndarray,
    lengths_um: np.ndarray,
    forward_mag: np.ndarray,
    reverse_mag: np.ndarray,
    width_um: float,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """Choose equal-width supports without looking at the geometry-edit outcome."""
    windows = candidate_windows(path_nodes, lengths_um, width_um)
    joint = np.asarray(forward_mag) * np.asarray(reverse_mag)

    def best_by(field):
        win = max(windows, key=lambda w: integrated_window_score(path_nodes, lengths_um, field, w))
        return support_from_window(path_nodes, win)

    out: dict[str, np.ndarray] = {
        "forward": best_by(forward_mag),
        "reverse": best_by(reverse_mag),
        "joint": best_by(joint),
    }
    random_win = windows[int(rng.integers(len(windows)))]
    out["random"] = support_from_window(path_nodes, random_win)

    joint_support = out["joint"]
    mids = path_edge_midpoints_um(path_nodes, lengths_um)
    idx_of = {int(node): i for i, node in enumerate(np.asarray(path_nodes, dtype=np.int64))}
    joint_center = float(
        np.average(
            [mids[idx_of[int(n)]] for n in joint_support],
            weights=np.asarray(lengths_um)[joint_support],
        )
    )

    # Distance-matched adversarial control: first require a grossly similar path
    # position and low overlap; if the joint support sits at an endpoint and no
    # distinct match exists, relax the position tolerance deterministically.
    chosen = None
    for position_tol in (0.25 * width_um, 0.5 * width_um, width_um, 2.0 * width_um):
        eligible = []
        for win in windows:
            s = support_from_window(path_nodes, win)
            center = float(
                np.average(
                    [mids[idx_of[int(n)]] for n in s],
                    weights=np.asarray(lengths_um)[s],
                )
            )
            if abs(center - joint_center) <= max(position_tol, 1.0):
                eligible.append((win, s, center))
        low_overlap = [row for row in eligible if overlap_fraction(row[1], joint_support, lengths_um) <= 0.25]
        pool = low_overlap if low_overlap else eligible
        if pool:
            chosen = min(
                pool,
                key=lambda row: integrated_window_score(path_nodes, lengths_um, joint, row[0]),
            )[1]
            if overlap_fraction(chosen, joint_support, lengths_um) < 0.999:
                break
    if chosen is None:
        chosen = out["random"].copy()
    out["distance_matched"] = chosen
    return out


def equal_budget_scales(support_cable_um: float, budget_um: float = EDIT_BUDGET_UM) -> tuple[float, float, float]:
    """Shorter/current/longer scales with equal absolute total cable change."""
    w = float(support_cable_um)
    if w <= 0:
        raise ValueError("support cable must be positive")
    frac = float(budget_um) / w
    if frac >= 0.95:
        raise ValueError("edit budget too large for support")
    return 1.0 - frac, 1.0, 1.0 + frac


def scaled_support_tones(tree, lengths, tones, tip, support, scale):
    trial = np.asarray(lengths, dtype=float).copy()
    trial[np.asarray(support, dtype=np.int64)] *= float(scale)
    return tone_transfer(tree, trial, np.asarray(tones, dtype=float), int(tip))


def bounded_score(z: np.ndarray, noise_power: float) -> float:
    p = np.abs(np.asarray(z, dtype=complex)) ** 2
    target = float(p[2])
    off = float(np.sum(p) - target)
    return target / max(off + float(noise_power), 1e-300)


def local_complex_electrotonic_increment(tree, lengths_um, omega_rad_s: float) -> np.ndarray:
    """Local cylindrical gamma*dx diagnostic for each morphology edge.

    gamma^2 = r_a_per_length * y_m_per_length.  The active admittance uses the
    exact same uniform quasi-active densities as Gate 13.  This is only a local
    surrogate coordinate for a tapered branched cable; the solved global operator
    remains the authority.
    """
    n = len(tree.parents)
    inc = np.zeros(n, dtype=complex)
    w = float(omega_rad_s)
    surface_admittance = (
        GLEAK_S_PER_CM2
        + 1j * w * CM_F_PER_CM2
        + GQ_S_PER_CM2 / (1.0 + 1j * w * TAU_Q_S)
    )
    for i in range(1, n):
        p = int(tree.parents[i])
        radius_um = max(0.5 * (float(tree.radii_um[p]) + float(tree.radii_um[i])), 0.05)
        a_cm = radius_um * 1.0e-4
        L_cm = max(float(lengths_um[i]), 1e-9) * 1.0e-4
        r_a = RA_OHM_CM / (np.pi * a_cm * a_cm)  # ohm / cm
        y_m = 2.0 * np.pi * a_cm * surface_admittance  # siemens / cm
        gamma = np.sqrt(r_a * y_m)
        if np.real(gamma) < 0:
            gamma = -gamma
        inc[i] = gamma * L_cm
    return inc


def evaluate_support(
    tree,
    lengths,
    tones,
    tip,
    support,
    z0,
    noise_power,
    omega_for_coordinate,
) -> dict:
    support = np.asarray(support, dtype=np.int64)
    cable = float(np.sum(np.asarray(lengths)[support]))
    scales = equal_budget_scales(cable)
    q0 = bounded_score(z0, noise_power)
    trial_rows = []
    best = None
    for scale in scales:
        if abs(scale - 1.0) < 1e-15:
            z = np.asarray(z0, dtype=complex)
        else:
            z = scaled_support_tones(tree, lengths, tones, tip, support, scale)
        q = bounded_score(z, noise_power)
        p = np.abs(z) ** 2
        scalar = np.vdot(z0, z) / max(float(np.vdot(z0, z0).real), 1e-300)
        residual = float(np.linalg.norm(z - scalar * z0) / max(np.linalg.norm(z), 1e-300))
        row = {
            "scale": float(scale),
            "signed_total_cable_change_um": float((scale - 1.0) * cable),
            "score": float(q),
            "delta_score": float(q - q0),
            "relative_delta_score": float((q - q0) / max(abs(q0), 1e-300)),
            "target_fraction": float(p[2] / max(float(np.sum(p)), 1e-300)),
            "target_gain_vs_baseline": float(abs(z[2]) / max(abs(z0[2]), 1e-300)),
            "target_phase_shift_rad": float(np.angle(z[2] / z0[2])),
            "five_tone_best_scalar_residual": residual,
        }
        trial_rows.append(row)
        if best is None or row["score"] > best["score"]:
            best = row
    assert best is not None
    gamma_inc = local_complex_electrotonic_increment(tree, lengths, omega_for_coordinate)
    electro_abs = float(np.sum(np.abs(gamma_inc[support])))
    electro_atten = float(np.sum(np.real(gamma_inc[support])))
    electro_phase = float(np.sum(np.abs(np.imag(gamma_inc[support]))))
    return {
        "support_nodes": int(len(support)),
        "support_first_node": int(support[0]),
        "support_last_node": int(support[-1]),
        "support_cable_um": cable,
        "equal_budget_total_abs_cable_change_um": float(EDIT_BUDGET_UM),
        "electrotonic_abs_gamma_width": electro_abs,
        "electrotonic_attenuation_width": electro_atten,
        "electrotonic_abs_phase_width": electro_phase,
        "trials": trial_rows,
        "gamma_delta_score": float(best["delta_score"]),
        "gamma_relative": float(best["relative_delta_score"]),
        "best_direction": "shorter" if best["scale"] < 1.0 else ("longer" if best["scale"] > 1.0 else "stay"),
        "best_scale": float(best["scale"]),
        "best_target_gain_vs_baseline": float(best["target_gain_vs_baseline"]),
        "best_non_scalar_residual": float(best["five_tone_best_scalar_residual"]),
    }


def carrier_setup(tree, lengths, tip, omega_target: float):
    tones = np.asarray([0.0, 0.35 * omega_target, omega_target, 2.5 * omega_target, 6.0 * omega_target])
    z0 = tone_transfer(tree, lengths, tones, tip)
    p0 = np.abs(z0) ** 2
    noise = float(OBSERVER_NOISE_FRACTION * p0[2])
    q0 = bounded_score(z0, noise)
    forward, ais = transfer_vector(tree, lengths, omega_target, tip)
    reverse, _ = transfer_vector(tree, lengths, omega_target, ais)
    return tones, z0, noise, q0, np.abs(forward[: len(tree.parents)]), np.abs(reverse[: len(tree.parents)]), ais


def run() -> dict:
    out_dir = Path("results")
    source = download_source(out_dir / "_source" / SOURCE_NAME)
    tree = load_dendritic_tree(source)
    lengths = edge_lengths_um(tree)
    tip, path, path_length = choose_distal_tip(tree, lengths)
    eligible_path = path[1:]  # each node names its parent->node edge

    # Define the baseline preferred carrier once, exactly as Gates 13/14 did.
    omegas = np.geomspace(0.5, 500.0, 120)
    gains = []
    for w in omegas:
        v, ais = transfer_vector(tree, lengths, float(w), tip)
        gains.append(abs(v[ais]))
    gains = np.asarray(gains, dtype=float)
    peak_i = int(np.argmax(gains))
    omega_star = float(omegas[peak_i])

    tones, z0, noise, q0, fmag, rmag, ais = carrier_setup(tree, lengths, tip, omega_star)
    rng = np.random.default_rng(RANDOM_SEED)

    selector_rows = []
    for width in WIDTHS_UM:
        supports = select_supports(eligible_path, lengths, fmag, rmag, float(width), rng)
        per_selector = {}
        for name in SELECTORS:
            per_selector[name] = evaluate_support(
                tree, lengths, tones, tip, supports[name], z0, noise, omega_star
            )
        joint = supports["joint"]
        per_selector["distance_matched"]["overlap_with_joint"] = overlap_fraction(
            supports["distance_matched"], joint, lengths
        )
        per_selector["random"]["overlap_with_joint"] = overlap_fraction(
            supports["random"], joint, lengths
        )
        selector_rows.append(
            {
                "requested_width_um": float(width),
                "selectors": per_selector,
            }
        )

    # Operator-native frequency attack: repeat only the joint selector at three
    # carrier frequencies, holding physical width and absolute edit budget fixed.
    frequency_rows = []
    for factor in FREQUENCY_FACTORS:
        omega = float(factor * omega_star)
        tones_f, zf, noise_f, qf, ff, rf, _ = carrier_setup(tree, lengths, tip, omega)
        rng_f = np.random.default_rng(RANDOM_SEED)  # irrelevant for joint, deterministic record
        rows = []
        for width in WIDTHS_UM:
            supports = select_supports(eligible_path, lengths, ff, rf, float(width), rng_f)
            ev = evaluate_support(tree, lengths, tones_f, tip, supports["joint"], zf, noise_f, omega)
            rows.append({"requested_width_um": float(width), **ev})
        gamma = np.asarray([row["gamma_relative"] for row in rows], dtype=float)
        max_gamma = float(np.max(gamma))
        if max_gamma > 0:
            half = 0.5 * max_gamma
            half_idx = int(np.flatnonzero(gamma >= half)[0])
            half_width_um = float(rows[half_idx]["support_cable_um"])
            half_electro = float(rows[half_idx]["electrotonic_abs_gamma_width"])
        else:
            half_width_um = None
            half_electro = None
        best_idx = int(np.argmax(gamma))
        frequency_rows.append(
            {
                "frequency_factor_vs_baseline_peak": float(factor),
                "omega_rad_s": omega,
                "baseline_bounded_score": float(qf),
                "rows": rows,
                "descriptive_half_max_onset_width_um": half_width_um,
                "descriptive_half_max_onset_electrotonic_abs_gamma": half_electro,
                "best_width_um": float(rows[best_idx]["support_cable_um"]),
                "best_gamma_relative": float(rows[best_idx]["gamma_relative"]),
            }
        )

    result = {
        "gate": 15,
        "question": (
            "With equal total cable-change budget, does operator leverage depend on the spatial support of a geometry edit, "
            "and is any apparent scale better described in micrometres or an operator-native electrotonic coordinate?"
        ),
        "morphology": {
            "dendritic_nodes_including_soma_root": int(len(tree.parents)),
            "selected_tip": int(tip),
            "selected_path_nodes": int(len(path)),
            "selected_path_length_um": float(path_length),
        },
        "baseline": {
            "omega_star_rad_s": omega_star,
            "peak_is_interior_to_frequency_sweep": bool(0 < peak_i < len(omegas) - 1),
            "observer_noise_fraction_of_baseline_target_power": float(OBSERVER_NOISE_FRACTION),
            "bounded_score": float(q0),
            "target_fraction": float(abs(z0[2]) ** 2 / max(float(np.sum(np.abs(z0) ** 2)), 1e-300)),
        },
        "protocol": {
            "requested_widths_um": [float(x) for x in WIDTHS_UM],
            "equal_budget_total_abs_cable_change_um": float(EDIT_BUDGET_UM),
            "selectors": list(SELECTORS),
            "frequency_factors": [float(x) for x in FREQUENCY_FACTORS],
            "gamma_definition": "best bounded-score improvement among shorter/current/longer edits with fixed sum |Delta L|",
            "distance_matched_definition": (
                "same gross path-position and physical width as joint where possible; choose low-joint, low-overlap contiguous support without using edit outcomes"
            ),
        },
        "selector_width_sweep_at_baseline_frequency": selector_rows,
        "joint_frequency_width_sweep": frequency_rows,
        "claim_boundary": (
            "This gate measures a support-response curve under one phenomenological quasi-active membrane and one real morphology. "
            "A rise, optimum, or half-max width is descriptive unless it survives resampling, branches, membrane parameters and other morphologies. "
            "The local electrotonic coordinate is a diagnostic derived from uniform-cylinder cable theory, not an assertion that the full branched active arbor has a literal wavelength."
        ),
    }
    out_dir.mkdir(exist_ok=True)
    (out_dir / "gate15_support_response.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    run()
