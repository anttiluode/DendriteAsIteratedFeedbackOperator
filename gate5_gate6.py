import json
from pathlib import Path

import numpy as np

from dendrite_operator import BranchedCable, effective_rank
from active_dendrite import (
    frequency_response,
    nmda_jacobian,
    resonant_state_matrix,
    simulate_nmda,
)


def run():
    cable = BranchedCable(branch_len=12)

    # Gate 5: a delayed restorative current turns the same physical branch
    # from a simple low-pass cable into a frequency-selective operator.
    passive_A = -(cable.G / cable.C[:, None])
    active_sites = cable.branch_b[cable.branch_len // 2 :]
    active_A = resonant_state_matrix(cable, active_sites, g_r=0.2, tau=20.0)

    omegas = np.concatenate([[0.0], np.logspace(-3.0, 0.5, 320)])
    tip_b = int(cable.branch_b[-1])
    passive_H = frequency_response(cable, passive_A, tip_b, cable.soma, omegas)
    active_H = frequency_response(cable, active_A, tip_b, cable.soma, omegas)
    passive_mag = np.abs(passive_H)
    active_mag = np.abs(active_H)

    passive_peak = int(np.argmax(passive_mag))
    active_peak = int(np.argmax(active_mag))
    poles = np.linalg.eigvals(active_A)
    complex_poles = poles[np.abs(poles.imag) > 1e-10]

    # Gate 6: local voltage-dependent conductances make the instantaneous
    # propagation operator depend on the current dendritic state.
    a_site = int(cable.branch_b[-3])
    b_site = int(cable.branch_b[-1])
    t0 = 5.0
    separation = 15.0
    nmda_kwargs = dict(
        g_peak=0.4,
        tau_s=8.0,
        vhalf=0.30,
        slope=0.08,
        reversal=1.0,
        dt=0.05,
        duration=60.0,
    )

    times, v_a = simulate_nmda(cable, [(a_site, t0)], **nmda_kwargs)
    _, v_b = simulate_nmda(cable, [(b_site, t0)], **nmda_kwargs)
    _, v_ab = simulate_nmda(cable, [(a_site, t0), (b_site, t0)], **nmda_kwargs)
    _, v_sep = simulate_nmda(
        cable, [(a_site, t0), (b_site, t0 + separation)], **nmda_kwargs
    )
    _, v_b_sep = simulate_nmda(cable, [(b_site, t0 + separation)], **nmda_kwargs)

    matched_residual = v_ab[:, cable.soma] - v_a[:, cable.soma] - v_b[:, cable.soma]
    separated_residual = (
        v_sep[:, cable.soma] - v_a[:, cable.soma] - v_b_sep[:, cable.soma]
    )
    matched_norm = float(np.linalg.norm(matched_residual))
    separated_norm = float(np.linalg.norm(separated_residual))

    # Hold the exact same external synaptic gates at t=8, but evaluate the
    # Jacobian once at the actual AB voltage state and once at V=0. This
    # isolates voltage-state dependence from a changed input schedule.
    eval_t = 8.0
    eval_idx = int(round(eval_t / nmda_kwargs["dt"]))
    events_ab = [(a_site, t0), (b_site, t0)]
    jac_kwargs = {
        k: nmda_kwargs[k]
        for k in ("g_peak", "tau_s", "vhalf", "slope", "reversal")
    }
    J_state = nmda_jacobian(cable, v_ab[eval_idx], eval_t, events_ab, **jac_kwargs)
    J_zero = nmda_jacobian(
        cable, np.zeros(cable.n), eval_t, events_ab, **jac_kwargs
    )
    delta_J = J_state - J_zero

    probe = np.zeros(cable.n)
    probe[b_site] = 1.0 / cable.C[b_site]
    readout = np.zeros(cable.n)
    readout[cable.soma] = 1.0
    h_zero = float(readout @ np.linalg.solve(-J_zero, probe))
    h_state = float(readout @ np.linalg.solve(-J_state, probe))

    max_real_eig = -np.inf
    for k in range(0, len(times), 20):
        if times[k] < t0:
            continue
        Jk = nmda_jacobian(cable, v_ab[k], times[k], events_ab, **jac_kwargs)
        max_real_eig = max(max_real_eig, float(np.max(np.linalg.eigvals(Jk).real)))

    receipt = {
        "model": {
            "compartments": cable.n,
            "branch_len": cable.branch_len,
        },
        "G5_active_frequency_selective_operator": {
            "active_sites": [int(x) for x in active_sites],
            "g_r": 0.2,
            "tau": 20.0,
            "passive_peak_omega": float(omegas[passive_peak]),
            "passive_peak_over_dc": float(passive_mag[passive_peak] / passive_mag[0]),
            "active_peak_omega": float(omegas[active_peak]),
            "active_peak_over_dc": float(active_mag[active_peak] / active_mag[0]),
            "active_max_pole_real_part": float(np.max(poles.real)),
            "number_of_complex_poles": int(len(complex_poles)),
            "largest_complex_pole_imag_abs": float(
                np.max(np.abs(complex_poles.imag)) if len(complex_poles) else 0.0
            ),
        },
        "G6_nmda_state_dependent_operator": {
            "a_site": a_site,
            "b_site": b_site,
            "matched_interaction_norm": matched_norm,
            "separated_interaction_norm": separated_norm,
            "matched_over_separated": float(matched_norm / separated_norm),
            "matched_peak_abs_residual": float(np.max(np.abs(matched_residual))),
            "separated_peak_abs_residual": float(np.max(np.abs(separated_residual))),
            "jacobian_eval_time": eval_t,
            "same_gate_state_vs_zero_deltaJ_norm": float(np.linalg.norm(delta_J)),
            "same_gate_state_vs_zero_deltaJ_rank": int(
                np.linalg.matrix_rank(delta_J, tol=1e-12)
            ),
            "instantaneous_B_probe_transfer_zero_state": h_zero,
            "instantaneous_B_probe_transfer_AB_state": h_state,
            "state_over_zero_probe_transfer": float(h_state / h_zero),
            "max_real_jacobian_eigenvalue_over_AB_run": max_real_eig,
            "locally_stable_in_this_run": bool(max_real_eig < 0.0),
        },
    }
    return receipt


if __name__ == "__main__":
    result = run()
    print(json.dumps(result, indent=2))
    out = Path("results")
    out.mkdir(exist_ok=True)
    (out / "gate5_gate6_receipt.json").write_text(json.dumps(result, indent=2) + "\n")
