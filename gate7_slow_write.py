import json
from pathlib import Path

import numpy as np

from dendrite_operator import BranchedCable, _add_edge
from active_dendrite import simulate_nmda


def step_operator_from_conductance(cable, G, dt=0.5):
    """Exact passive step for an edited conductance matrix."""
    c_isqrt = 1.0 / np.sqrt(cable.C)
    B = (c_isqrt[:, None] * G) * c_isqrt[None, :]
    lam, Q = np.linalg.eigh(B)
    E = (Q * np.exp(-float(dt) * lam)) @ Q.T
    c_sqrt = np.sqrt(cable.C)
    return (c_isqrt[:, None] * E) * c_sqrt[None, :]


def probe_trace(cable, G, site, steps=180, dt=0.5):
    """Erase all fast state, inject one unit charge, and follow the passive response."""
    P = step_operator_from_conductance(cable, G, dt=dt)
    x = np.zeros(cable.n)
    x[int(site)] = 1.0 / cable.C[int(site)]
    history = np.empty((steps, cable.n))
    for k in range(steps):
        history[k] = x
        x = P @ x
    return history


def local_slow_write(voltage_history, write_site, threshold=0.012, eta=800.0, dt=0.05):
    """Local activity-dependent conductance write.

    The rule intentionally has no knowledge of A or B identity:
        Delta g = eta * integral max(V_write - threshold, 0)^2 dt.

    It is a toy plasticity rule, not a fitted biological learning law.
    """
    excess = np.maximum(voltage_history[:, int(write_site)] - float(threshold), 0.0)
    return float(eta * np.sum(excess * excess) * float(dt))


def run():
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

    _, v0 = simulate_nmda(cable, [], **sim_kwargs)
    _, va = simulate_nmda(cable, [(a_site, t0)], **sim_kwargs)
    _, vb = simulate_nmda(cable, [(b_site, t0)], **sim_kwargs)
    _, vab = simulate_nmda(cable, [(a_site, t0), (b_site, t0)], **sim_kwargs)
    _, vsep = simulate_nmda(
        cable, [(a_site, t0), (b_site, t0 + separation)], **sim_kwargs
    )

    write_kwargs = dict(
        write_site=write_site,
        threshold=0.012,
        eta=800.0,
        dt=sim_kwargs["dt"],
    )
    writes = {
        "W0": local_slow_write(v0, **write_kwargs),
        "WA": local_slow_write(va, **write_kwargs),
        "WB": local_slow_write(vb, **write_kwargs),
        "WAB": local_slow_write(vab, **write_kwargs),
        "Wsep": local_slow_write(vsep, **write_kwargs),
    }
    isolated = writes["WAB"] - writes["WA"] - writes["WB"] + writes["W0"]

    # The slow variable edits one local axial constraint. Every fast voltage
    # and synaptic gate is then discarded: recall starts from x=0.
    edge_i = write_site
    edge_j = b_site

    traces = {}
    for world, delta_g in writes.items():
        G = cable.G.copy()
        _add_edge(G, edge_i, edge_j, delta_g)
        traces[world] = probe_trace(cable, G, a_site, steps=180, dt=0.5)

    base = traces["W0"]
    ab = traces["WAB"]
    sep = traces["Wsep"]

    base_b_peak = float(np.max(base[:, b_site]))
    ab_b_peak = float(np.max(ab[:, b_site]))
    sep_b_peak = float(np.max(sep[:, b_site]))
    base_soma_peak = float(np.max(base[:, cable.soma]))
    ab_soma_peak = float(np.max(ab[:, cable.soma]))
    sep_soma_peak = float(np.max(sep[:, cable.soma]))

    rel_b = (ab_b_peak - base_b_peak) / base_b_peak
    rel_soma = (ab_soma_peak - base_soma_peak) / base_soma_peak
    rel_sep_b = (sep_b_peak - base_b_peak) / base_b_peak

    receipt = {
        "model": {
            "compartments": cable.n,
            "a_site": a_site,
            "write_site": write_site,
            "b_site": b_site,
            "edited_edge": [edge_i, edge_j],
        },
        "G7_persistent_local_write": {
            "write_rule": "eta * integral(max(V_write-threshold,0)^2) dt",
            "threshold": write_kwargs["threshold"],
            "eta": write_kwargs["eta"],
            "delta_g_by_world": writes,
            "collision_specific_delta_g": isolated,
            "matched_over_separated_write": (
                float(writes["WAB"] / writes["Wsep"])
                if writes["Wsep"] > 0
                else float("inf")
            ),
            "individual_A_write": writes["WA"],
            "individual_B_write": writes["WB"],
        },
        "G7_after_fast_state_erasure": {
            "recall_initial_fast_state": "exact zero",
            "recall_synaptic_gates": "none; passive unit-charge A probe only",
            "B_peak_baseline": base_b_peak,
            "B_peak_matched_memory": ab_b_peak,
            "B_peak_separated_memory": sep_b_peak,
            "relative_B_change_matched": float(rel_b),
            "relative_B_change_separated": float(rel_sep_b),
            "soma_peak_baseline": base_soma_peak,
            "soma_peak_matched_memory": ab_soma_peak,
            "soma_peak_separated_memory": sep_soma_peak,
            "relative_soma_change_matched": float(rel_soma),
            "branch_vs_soma_relative_visibility": float(abs(rel_b / rel_soma)),
        },
    }
    return receipt


if __name__ == "__main__":
    result = run()
    print(json.dumps(result, indent=2))
    out = Path("results")
    out.mkdir(exist_ok=True)
    (out / "slow_write_receipt.json").write_text(json.dumps(result, indent=2) + "\n")
