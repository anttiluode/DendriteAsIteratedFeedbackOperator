"""Diagnostic: where along branch A does Gate-5 mode purification exist locally?

No nonlinear threshold is chosen here.  We simply measure the frequency-domain
transfer from the distal tip to every branch compartment for the active and
matched-passive substrates.  This prevents the next nonlinear gate from
quietly choosing a hotspot before we know what the distributed operator itself
does there.
"""

from __future__ import annotations

import json

import numpy as np

from active_resonance import QuasiActiveDendrite, log_frequency_grid, peak_info


def target_fraction(gains, target_index=2):
    p = np.square(np.asarray(gains, dtype=float))
    return float(p[target_index] / p.sum())


def run():
    model = QuasiActiveDendrite()
    tip = int(model.branch_a[-1])
    ws = log_frequency_grid(0.002, 1.0, 600)

    rows = []
    for site in model.branch_a:
        site = int(site)
        spectrum = model.spectrum(ws, tip, site)
        peak = peak_info(ws, spectrum)
        w = peak["omega"]
        # Exclude DC here.  We want a true multimode temporal mixture rather
        # than allowing the passive cable's large DC response to settle the
        # question by itself.
        tones = np.asarray([0.35 * w, 0.60 * w, w, 2.5 * w, 6.0 * w])
        active = model.spectrum(tones, tip, site)
        passive = model.passive_spectrum(tones, tip, site)
        rows.append({
            "site": site,
            "distance_from_soma_index": int(np.where(model.branch_a == site)[0][0] + 1),
            "active_peak_omega": float(w),
            "active_peak_over_dc": float(peak["peak_over_dc"]),
            "active_target_fraction": target_fraction(active),
            "passive_target_fraction_same_tones": target_fraction(passive),
            "active_over_passive_target_fraction": float(target_fraction(active) / max(target_fraction(passive), 1e-15)),
            "active_total_gain_power": float(np.sum(active**2)),
            "passive_total_gain_power": float(np.sum(passive**2)),
            "tones": [float(x) for x in tones],
            "active_gains": [float(x) for x in active],
            "passive_gains": [float(x) for x in passive],
        })

    best = max(rows, key=lambda r: r["active_target_fraction"])
    best_contrast = max(rows, key=lambda r: r["active_over_passive_target_fraction"])
    return {
        "tip": tip,
        "rows": rows,
        "best_active_purity": best,
        "best_active_vs_passive_purity_contrast": best_contrast,
        "note": "This is a diagnostic only. No nonlinear threshold or preferred hotspot is fitted here."
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
