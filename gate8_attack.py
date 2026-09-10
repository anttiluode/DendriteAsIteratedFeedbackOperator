import json
from pathlib import Path

import numpy as np

from dendrite_operator import BranchedCable, _add_edge
from active_dendrite import simulate_nmda
from gate7_slow_write import probe_trace


def _raw_write(history, site, threshold, dt):
    excess = np.maximum(history[:, int(site)] - float(threshold), 0.0)
    return float(np.sum(excess * excess) * float(dt))


def _route_change(cable, a_site, b_site, edge_i, edge_j, delta_g, baseline_peak):
    G = cable.G.copy()
    _add_edge(G, edge_i, edge_j, float(delta_g))
    trace = probe_trace(cable, G, a_site, steps=180, dt=0.5)
    peak = float(np.max(trace[:, b_site]))
    return float((peak - baseline_peak) / baseline_peak)


def run():
    cable = BranchedCable(branch_len=12)
    a_site = int(cable.branch_b[-3])
    write_site = int(cable.branch_b[-2])
    b_site = int(cable.branch_b[-1])
    t0 = 5.0

    sim_kwargs = dict(
        g_peak=0.4,
        tau_s=8.0,
        vhalf=0.30,
        slope=0.08,
        reversal=1.0,
        dt=0.05,
        duration=60.0,
    )

    _, va = simulate_nmda(cable, [(a_site, t0)], **sim_kwargs)
    _, vb = simulate_nmda(cable, [(b_site, t0)], **sim_kwargs)
    _, vab = simulate_nmda(cable, [(a_site, t0), (b_site, t0)], **sim_kwargs)
    _, vsep = simulate_nmda(
        cable, [(a_site, t0), (b_site, t0 + 15.0)], **sim_kwargs
    )

    baseline = probe_trace(cable, cable.G.copy(), a_site, steps=180, dt=0.5)
    baseline_b_peak = float(np.max(baseline[:, b_site]))

    eta = 800.0
    thresholds = np.arange(0.007, 0.0155 + 1e-12, 0.0005)
    threshold_rows = []

    for threshold in thresholds:
        raw_a = _raw_write(va, write_site, threshold, sim_kwargs["dt"])
        raw_b = _raw_write(vb, write_site, threshold, sim_kwargs["dt"])
        raw_ab = _raw_write(vab, write_site, threshold, sim_kwargs["dt"])
        raw_sep = _raw_write(vsep, write_site, threshold, sim_kwargs["dt"])

        wa = eta * raw_a
        wb = eta * raw_b
        wab = eta * raw_ab
        wsep = eta * raw_sep

        contamination = (wa + wb) / wab if wab > 0 else float("inf")
        write_ratio = wab / wsep if wsep > 0 else float("inf")
        route_ab = _route_change(
            cable, a_site, b_site, write_site, b_site, wab, baseline_b_peak
        )
        route_sep = _route_change(
            cable, a_site, b_site, write_site, b_site, wsep, baseline_b_peak
        )

        robust = bool(
            wab > 0
            and contamination <= 0.01
            and write_ratio >= 3.0
            and route_ab >= 0.05
            and route_sep <= 0.02
        )
        threshold_rows.append(
            {
                "threshold": float(threshold),
                "WA": float(wa),
                "WB": float(wb),
                "WAB": float(wab),
                "Wsep": float(wsep),
                "individual_contamination_fraction": float(contamination),
                "matched_over_separated_write": float(write_ratio),
                "matched_route_change": float(route_ab),
                "separated_route_change": float(route_sep),
                "robust": robust,
            }
        )

    robust_thresholds = [r["threshold"] for r in threshold_rows if r["robust"]]

    fixed_threshold = 0.012
    raw_ab = _raw_write(vab, write_site, fixed_threshold, sim_kwargs["dt"])
    raw_sep = _raw_write(vsep, write_site, fixed_threshold, sim_kwargs["dt"])
    eta_values = [100.0, 200.0, 400.0, 800.0, 1200.0, 1600.0]
    eta_rows = []
    for eta_i in eta_values:
        wab = eta_i * raw_ab
        wsep = eta_i * raw_sep
        route_ab = _route_change(
            cable, a_site, b_site, write_site, b_site, wab, baseline_b_peak
        )
        route_sep = _route_change(
            cable, a_site, b_site, write_site, b_site, wsep, baseline_b_peak
        )
        route_ratio = route_ab / route_sep if route_sep > 0 else float("inf")
        eta_rows.append(
            {
                "eta": eta_i,
                "WAB": float(wab),
                "Wsep": float(wsep),
                "matched_route_change": float(route_ab),
                "separated_route_change": float(route_sep),
                "matched_over_separated_route_change": float(route_ratio),
            }
        )

    receipt = {
        "G8_threshold_attack": {
            "thresholds_tested": int(len(threshold_rows)),
            "criterion": (
                "individual contamination <=1%, matched/separated write >=3x, "
                "matched B-route change >=5%, separated route change <=2%"
            ),
            "robust_threshold_count": int(len(robust_thresholds)),
            "robust_threshold_min": (
                float(min(robust_thresholds)) if robust_thresholds else None
            ),
            "robust_threshold_max": (
                float(max(robust_thresholds)) if robust_thresholds else None
            ),
            "frozen_gate7_threshold": fixed_threshold,
            "frozen_threshold_is_in_robust_window": bool(
                robust_thresholds
                and min(robust_thresholds) - 1e-12 <= fixed_threshold
                <= max(robust_thresholds) + 1e-12
            ),
            "rows": threshold_rows,
        },
        "G8_scale_attack": {
            "fixed_threshold": fixed_threshold,
            "eta_values": eta_values,
            "minimum_route_selectivity_ratio": float(
                min(r["matched_over_separated_route_change"] for r in eta_rows)
            ),
            "maximum_route_selectivity_ratio": float(
                max(r["matched_over_separated_route_change"] for r in eta_rows)
            ),
            "rows": eta_rows,
        },
    }
    return receipt


if __name__ == "__main__":
    result = run()
    print(json.dumps(result, indent=2))
    out = Path("results")
    out.mkdir(exist_ok=True)
    (out / "gate8_attack_receipt.json").write_text(json.dumps(result, indent=2) + "\n")
