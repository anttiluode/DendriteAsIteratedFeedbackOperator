import json
from pathlib import Path

import numpy as np

from dendrite_operator import BranchedCable, _add_edge
from active_dendrite import simulate_nmda
from gate7_slow_write import probe_trace


def energy_write(history, site, eta, dt):
    """Threshold-free local slow write: eta * integral V(site,t)^2 dt."""
    v = history[:, int(site)]
    return float(eta * np.sum(v * v) * float(dt))


def run(eta=30.0):
    cable = BranchedCable(branch_len=12)
    a_site = int(cable.branch_b[-3])
    write_site = int(cable.branch_b[-2])
    b_site = int(cable.branch_b[-1])
    t0 = 5.0
    separation = 15.0

    sim_kwargs = dict(
        g_peak=0.4,
        tau_s=8.0,
        vhalf=0.30,
        slope=0.08,
        reversal=1.0,
        dt=0.05,
        duration=60.0,
    )

    events = {
        "W0": [],
        "WA": [(a_site, t0)],
        "WB": [(b_site, t0)],
        "WAB": [(a_site, t0), (b_site, t0)],
        "Wsep": [(a_site, t0), (b_site, t0 + separation)],
        "WBsep": [(b_site, t0 + separation)],
    }

    voltage = {}
    for world, ev in events.items():
        _, voltage[world] = simulate_nmda(cable, ev, **sim_kwargs)

    writes = {
        world: energy_write(v, write_site, eta=eta, dt=sim_kwargs["dt"])
        for world, v in voltage.items()
    }

    isolated_write_matched = (
        writes["WAB"] - writes["WA"] - writes["WB"] + writes["W0"]
    )
    isolated_write_separated = (
        writes["Wsep"] - writes["WA"] - writes["WBsep"] + writes["W0"]
    )

    # Each world keeps its own total physical write.  We do NOT install the
    # isolated counterfactual component into a world.  After exact fast-state
    # erasure, every world receives the same passive A probe.  Inclusion/
    # exclusion is applied only to the measured later responses.
    traces = {}
    for world, delta_g in writes.items():
        G = cable.G.copy()
        _add_edge(G, write_site, b_site, delta_g)
        traces[world] = probe_trace(cable, G, a_site, steps=180, dt=0.5)

    chi_matched = (
        traces["WAB"] - traces["WA"] - traces["WB"] + traces["W0"]
    )
    chi_separated = (
        traces["Wsep"] - traces["WA"] - traces["WBsep"] + traces["W0"]
    )

    chi_b = chi_matched[:, b_site]
    chi_sep_b = chi_separated[:, b_site]
    chi_soma = chi_matched[:, cable.soma]
    chi_sep_soma = chi_separated[:, cable.soma]

    # Attack the arbitrary scale: inclusion/exclusion selectivity should not
    # disappear merely because eta changes.
    eta_sweep = [10.0, 20.0, 30.0, 50.0, 100.0]
    sweep_rows = []
    raw_writes = {
        world: energy_write(v, write_site, eta=1.0, dt=sim_kwargs["dt"])
        for world, v in voltage.items()
    }
    for eta_i in eta_sweep:
        world_traces = {}
        for world, raw in raw_writes.items():
            G = cable.G.copy()
            _add_edge(G, write_site, b_site, eta_i * raw)
            world_traces[world] = probe_trace(cable, G, a_site, steps=180, dt=0.5)
        cm = (
            world_traces["WAB"]
            - world_traces["WA"]
            - world_traces["WB"]
            + world_traces["W0"]
        )[:, b_site]
        cs = (
            world_traces["Wsep"]
            - world_traces["WA"]
            - world_traces["WBsep"]
            + world_traces["W0"]
        )[:, b_site]
        sweep_rows.append(
            {
                "eta": eta_i,
                "matched_chi_norm": float(np.linalg.norm(cm)),
                "separated_chi_norm": float(np.linalg.norm(cs)),
                "matched_over_separated_chi_norm": float(
                    np.linalg.norm(cm) / np.linalg.norm(cs)
                ),
            }
        )

    receipt = {
        "model": {
            "compartments": cable.n,
            "a_site": a_site,
            "write_site": write_site,
            "b_site": b_site,
            "separation": separation,
        },
        "G9_threshold_free_write": {
            "rule": "eta * integral(V_write^2) dt",
            "eta": eta,
            "delta_g_by_world": writes,
            "individual_A_write": writes["WA"],
            "individual_B_write": writes["WB"],
            "collision_specific_write_matched": float(isolated_write_matched),
            "collision_specific_write_separated": float(isolated_write_separated),
            "matched_over_separated_collision_write": float(
                isolated_write_matched / isolated_write_separated
            ),
        },
        "G9_after_fast_state_erasure": {
            "recall_initial_fast_state": "exact zero in every world",
            "recall_probe": "same passive unit-charge A cue in every world",
            "matched_B_chi_norm": float(np.linalg.norm(chi_b)),
            "separated_B_chi_norm": float(np.linalg.norm(chi_sep_b)),
            "matched_over_separated_B_chi_norm": float(
                np.linalg.norm(chi_b) / np.linalg.norm(chi_sep_b)
            ),
            "matched_B_chi_peak_abs": float(np.max(np.abs(chi_b))),
            "separated_B_chi_peak_abs": float(np.max(np.abs(chi_sep_b))),
            "matched_over_separated_B_chi_peak": float(
                np.max(np.abs(chi_b)) / np.max(np.abs(chi_sep_b))
            ),
            "matched_soma_chi_norm": float(np.linalg.norm(chi_soma)),
            "separated_soma_chi_norm": float(np.linalg.norm(chi_sep_soma)),
        },
        "G9_eta_attack": {
            "eta_values": eta_sweep,
            "minimum_B_chi_selectivity": float(
                min(r["matched_over_separated_chi_norm"] for r in sweep_rows)
            ),
            "maximum_B_chi_selectivity": float(
                max(r["matched_over_separated_chi_norm"] for r in sweep_rows)
            ),
            "rows": sweep_rows,
        },
        "interpretation": (
            "Inclusion/exclusion is an analysis across counterfactual worlds, not a "
            "mechanism available to one neuron. It tests whether a persistent later "
            "response contains a component not explained by adding the two individual writes."
        ),
    }
    return receipt


if __name__ == "__main__":
    result = run()
    print(json.dumps(result, indent=2, allow_nan=False))
    out = Path("results")
    out.mkdir(exist_ok=True)
    (out / "gate9_causal_write_receipt.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
