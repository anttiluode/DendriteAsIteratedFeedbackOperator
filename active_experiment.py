"""Gate 5: does an active dendritic operator actually purify a non-zero mode?

The passive gates established only diffusive mode selection.  This gate asks a
harder question: can a dendrite-like cable with a restorative active membrane
have a non-zero pass band, complex poles, and enrich a selected temporal mode
from an equal-tone mixture?
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from active_resonance import (
    QuasiActiveDendrite,
    equal_tone_purification,
    log_frequency_grid,
    peak_info,
)
from dendrite_operator import effective_rank


def _pole_summary(poles):
    poles = np.asarray(poles)
    complex_mask = np.abs(poles.imag) > 1e-8
    cp = poles[complex_mask]
    least_damped = None
    if len(cp):
        p = cp[int(np.argmax(cp.real))]
        least_damped = {
            "real": float(p.real),
            "imag_abs": float(abs(p.imag)),
            "decay_time": float(-1.0 / p.real) if p.real < 0 else None,
        }
    return {
        "max_real_part": float(np.max(poles.real)),
        "stable": bool(np.max(poles.real) < 0.0),
        "complex_pole_count": int(np.sum(complex_mask)),
        "least_damped_complex_pole": least_damped,
    }


def _rank_one_edit_receipt(model, omega, site, delta_g=0.06):
    H0 = model.transfer_matrix(omega)
    H1 = model.edited_transfer_matrix(omega, site, delta_g)
    dH = H1 - H0
    s = np.linalg.svd(dH, compute_uv=False)
    return {
        "site": int(site),
        "omega": float(omega),
        "delta_g": float(delta_g),
        "deltaH_effective_rank": effective_rank(s),
        "top_singular_fraction": float(s[0] / s.sum()),
        "dense_fraction_above_2pct_peak": float(np.mean(np.abs(dH) > 0.02 * np.max(np.abs(dH)))),
        "s1": float(s[0]),
        "s2_over_s1": float(s[1] / s[0]) if len(s) > 1 else 0.0,
    }


def run():
    model = QuasiActiveDendrite(branch_len=12)
    omegas = log_frequency_grid(0.002, 1.0, 600)

    tip_a = int(model.branch_a[-1])
    tip_b = int(model.branch_b[-1])

    active_a = model.spectrum(omegas, tip_a, model.soma)
    active_b = model.spectrum(omegas, tip_b, model.soma)
    passive_a = model.passive_spectrum(omegas, tip_a, model.soma)
    passive_b = model.passive_spectrum(omegas, tip_b, model.soma)

    pa = peak_info(omegas, active_a)
    pb = peak_info(omegas, active_b)
    ppa = peak_info(omegas, passive_a)
    ppb = peak_info(omegas, passive_b)

    pur_a = equal_tone_purification(model, tip_a, pa["omega"], model.soma)
    pur_b = equal_tone_purification(model, tip_b, pb["omega"], model.soma)

    pole = _pole_summary(model.poles())

    # A local active-membrane edit remains a local constraint edit.  At a fixed
    # temporal frequency it is a rank-one diagonal admittance perturbation and
    # therefore produces a dense rank-one resolvent change globally.
    edit_site = int(model.branch_a[len(model.branch_a) // 2])
    edit = _rank_one_edit_receipt(model, pa["omega"], edit_site)

    # Functional tuning consequence at the same distal-to-soma channel.
    H_before = model.transfer_matrix(pa["omega"])
    H_after = model.edited_transfer_matrix(pa["omega"], edit_site, 0.06)
    gain_before = abs(H_before[model.soma, tip_a])
    gain_after = abs(H_after[model.soma, tip_a])
    edit["tipA_to_soma_gain_before"] = float(gain_before)
    edit["tipA_to_soma_gain_after"] = float(gain_after)
    edit["relative_gain_change"] = float((gain_after - gain_before) / gain_before)

    receipt = {
        "model": {
            "compartments": model.n,
            "quasi_active_state_variables": int(len(model.active_sites)),
            "branch_A": {"g_q": 0.30, "tau_q": 18.0},
            "branch_B": {"g_q": 0.18, "tau_q": 42.0},
            "interpretation": "phenomenological restorative quasi-active membrane; not a fitted HCN model",
        },
        "G5a_passive_control": {
            "branch_A_peak": ppa,
            "branch_B_peak": ppb,
            "expectation": "passive control should peak at DC",
        },
        "G5b_active_nonzero_resonance": {
            "branch_A_peak": pa,
            "branch_B_peak": pb,
            "peak_frequency_ratio_A_over_B": float(pa["omega"] / pb["omega"]) if pb["omega"] > 0 else None,
            "both_peaks_nonzero": bool(pa["omega"] > 0 and pb["omega"] > 0),
        },
        "G5c_complex_poles": pole,
        "G5d_equal_tone_mode_purification": {
            "branch_A": pur_a,
            "branch_B": pur_b,
            "definition": "target temporal mode carries 20% of equal-tone input energy; purification means a larger target fraction at the bounded output",
        },
        "G5e_local_active_edit_global_operator_change": edit,
    }

    # Gate verdicts are relationship tests, not hard-coded desired numbers.
    receipt["verdict"] = {
        "passive_is_dc_lowpass": bool(ppa["omega"] == 0.0 and ppb["omega"] == 0.0),
        "active_has_nonzero_passbands": bool(pa["omega"] > 0.0 and pb["omega"] > 0.0),
        "active_has_complex_stable_poles": bool(pole["stable"] and pole["complex_pole_count"] >= 2),
        "branch_A_enriches_target": bool(pur_a["active_output_target_energy_fraction"] > 0.2),
        "branch_B_enriches_target": bool(pur_b["active_output_target_energy_fraction"] > 0.2),
        "local_edit_is_rank_one": bool(edit["s2_over_s1"] < 1e-9),
    }
    receipt["verdict"]["passes"] = bool(all(receipt["verdict"].values()))
    return receipt


if __name__ == "__main__":
    r = run()
    print(json.dumps(r, indent=2))
    out = Path("results")
    out.mkdir(exist_ok=True)
    (out / "active_receipt.json").write_text(json.dumps(r, indent=2) + "\n")
    if not r["verdict"]["passes"]:
        raise SystemExit("Gate 5 relationship verdict failed")
