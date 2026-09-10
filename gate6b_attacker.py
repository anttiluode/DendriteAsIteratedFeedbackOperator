"""Gate 6b attacker: is Gate 6 just a hand-placed threshold and phasor trick?

Two attacks are deliberately *not* written as pass/fail tests.

1. Remove all quasi-active membrane conductance, but keep the exact same raw
   sources, carrier, phases, amplitudes, hotspot, nonlinear current, and knee.
   This asks whether the resonant substrate itself matters.
2. Sweep the nonlinear knee without recalibrating anything. This asks whether
   the matched advantage exists only at the one threshold we chose.

The script prints and freezes the evidence. An embarrassing result is allowed.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from active_resonance import QuasiActiveDendrite
from nonlinear_sieve import ResonantNonlinearSieve


def rms(x):
    return float(np.sqrt(np.mean(np.square(x))))


def case_residual(machine, spec):
    soma_lin, hot_lin = machine.simulate(spec, nonlinear=False)
    soma_nl, _ = machine.simulate(spec, nonlinear=True)
    return {
        "linear_hotspot_peak": float(np.max(np.abs(hot_lin))),
        "linear_soma_rms": rms(soma_lin),
        "nonlinear_residual_soma_rms": rms(soma_nl - soma_lin),
    }


def threshold_sweep(machine, specs, thresholds):
    # Linear references do not depend on the nonlinear knee, so compute once.
    refs = {}
    for name in ("matched", "quadrature"):
        soma_lin, hot_lin = machine.simulate(specs[name], nonlinear=False)
        refs[name] = (soma_lin, hot_lin)

    original = machine.threshold
    rows = []
    try:
        for threshold in thresholds:
            machine.threshold = float(threshold)
            vals = {}
            for name in ("matched", "quadrature"):
                soma_lin, hot_lin = refs[name]
                soma_nl, _ = machine.simulate(specs[name], nonlinear=True)
                vals[name] = {
                    "linear_hotspot_peak": float(np.max(np.abs(hot_lin))),
                    "nonlinear_residual_soma_rms": rms(soma_nl - soma_lin),
                }
            m = vals["matched"]["nonlinear_residual_soma_rms"]
            q = vals["quadrature"]["nonlinear_residual_soma_rms"]
            rows.append({
                "threshold": float(threshold),
                "matched": vals["matched"],
                "quadrature": vals["quadrature"],
                "matched_over_quadrature": float(m / max(q, 1e-15)),
            })
    finally:
        machine.threshold = original
    return rows


def run():
    active = ResonantNonlinearSieve()
    specs = active.case_specs()

    # Exact protocol ablation: only the quasi-active conductances disappear.
    passive_model = QuasiActiveDendrite(g_a=0.0, g_b=0.0)
    passive = active.same_protocol_on(passive_model)

    names = ("matched", "quadrature", "wrong_high_carrier", "spatially_separated", "single_A")
    active_cases = {name: case_residual(active, specs[name]) for name in names}
    passive_cases = {name: case_residual(passive, specs[name]) for name in names}

    def ratio(cases, a, b):
        x = cases[a]["nonlinear_residual_soma_rms"]
        y = cases[b]["nonlinear_residual_soma_rms"]
        return float(x / max(y, 1e-15))

    thresholds = [0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075]
    sweep = threshold_sweep(active, specs, thresholds)
    ratios = [row["matched_over_quadrature"] for row in sweep]

    return {
        "question": "Does Gate 6 survive removal of quasi-active resonance and movement of the nonlinear knee?",
        "fixed_protocol": {
            "omega": active.omega,
            "phase_b_matched": active.phase_b_matched,
            "relative_input_amplitudes": [active.rel_amp_a, active.rel_amp_b],
            "drive_scale": active.drive_scale,
            "hotspot": active.hotspot,
            "default_threshold": active.threshold,
        },
        "active_substrate": active_cases,
        "passive_ablation_same_protocol": passive_cases,
        "comparisons": {
            "active_matched_over_quadrature": ratio(active_cases, "matched", "quadrature"),
            "passive_matched_over_quadrature": ratio(passive_cases, "matched", "quadrature"),
            "active_matched_over_wrong_high": ratio(active_cases, "matched", "wrong_high_carrier"),
            "passive_matched_over_wrong_high": ratio(passive_cases, "matched", "wrong_high_carrier"),
            "active_matched_residual": active_cases["matched"]["nonlinear_residual_soma_rms"],
            "passive_matched_residual": passive_cases["matched"]["nonlinear_residual_soma_rms"],
            "active_to_passive_matched_residual": float(
                active_cases["matched"]["nonlinear_residual_soma_rms"]
                / max(passive_cases["matched"]["nonlinear_residual_soma_rms"], 1e-15)
            ),
            "active_to_passive_matched_hotspot_peak": float(
                active_cases["matched"]["linear_hotspot_peak"]
                / max(passive_cases["matched"]["linear_hotspot_peak"], 1e-15)
            ),
        },
        "active_threshold_sweep": sweep,
        "threshold_sweep_summary": {
            "min_matched_over_quadrature": float(min(ratios)),
            "max_matched_over_quadrature": float(max(ratios)),
            "all_ratios_above_one": bool(all(r > 1.0 for r in ratios)),
            "note": "This sweep changes only the nonlinear knee; carrier, phases, amplitudes and substrate stay frozen."
        },
        "interpretation_rule": {
            "resonance_is_not_necessary_if": "passive ablation retains comparable nonlinear selectivity under the identical raw protocol",
            "single_threshold_story_is_weakened_if": "matched/quadrature advantage persists across the knee sweep",
            "warning": "Even a positive result remains a designed mechanism proof, not evidence that a biological dendrite discovered these parameters."
        }
    }


if __name__ == "__main__":
    receipt = run()
    print(json.dumps(receipt, indent=2))
    out = Path("results")
    out.mkdir(exist_ok=True)
    (out / "gate6b_attacker.json").write_text(json.dumps(receipt, indent=2) + "\n")
