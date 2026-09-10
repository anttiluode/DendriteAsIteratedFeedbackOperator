"""Gate 12: a geometric length operator, plus the silence-trap attack.

The first run exposed the important failure: maximizing *purity alone* makes the
chosen dendritic segment grow without bound.  The target fraction rises while
the absolute transmitted signal disappears.  That is the same conceptual trap
as an ungrounded SighImageSuper mode: a direction can dominate because
Everything Else died, not because a useful answer remains.

Gate 12 therefore separates two questions.

1. PURITY-ONLY length update
       L -> whichever nearby L gives the largest target fraction.
   This has no finite fixed point in the current toy.

2. GROUNDED length update for a bounded observer
       score = target_power / (off_target_power + observer_noise_floor).
   Now an infinitely isolated branch scores zero because the observer cannot
   recover a vanishing target through its own noise floor.  A finite length can
   become a fixed point of the local geometry update.

The geometry proxy remains deliberately simple.  For one compartment segment,
length L -> s L at fixed radius and membrane/channel density implies

    axial conductance      g_ax -> g_ax / s
    membrane capacitance  C    -> s C
    leak conductance      g_L  -> s g_L
    quasi-active membrane admittance -> s times local value

This is compartmental cable physics, not a literal acoustic wavelength model.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from active_resonance import log_frequency_grid
from gate11_ais_receipt import DendriteAIS, _edge_outer


class LengthOperator:
    """One-segment length edit applied to the Gate-11 dendrite+AIS operator."""

    def __init__(self, machine: DendriteAIS, parent: int, site: int):
        self.machine = machine
        self.parent = int(parent)
        self.site = int(site)
        if self.site >= machine.n_d or self.parent >= machine.n_d:
            raise ValueError("Gate 12 currently edits a dendritic segment only")

        self.g_ax0 = float(-machine.dendrite.G[self.parent, self.site])
        if self.g_ax0 <= 0:
            raise ValueError("selected parent/site are not an axial edge")
        self.c0 = float(machine.dendrite.C[self.site])
        self.g_leak0 = float(machine.dendrite.base.g_leak[self.site])

    def admittance(self, omega: float, length_scale: float) -> np.ndarray:
        s = float(length_scale)
        if s <= 0:
            raise ValueError("length_scale must be positive")
        w = float(omega)
        Y = self.machine.admittance(w).copy()

        # Replace original axial edge g by g/s.
        Y += self.g_ax0 * (1.0 / s - 1.0) * _edge_outer(
            self.machine.n, self.parent, self.site
        )

        # At fixed radius/density, local membrane area scales with segment length.
        q0 = complex(self.machine.dendrite.quasi_admittance(w)[self.site])
        membrane0 = self.g_leak0 + 1j * w * self.c0 + q0
        Y[self.site, self.site] += (s - 1.0) * membrane0
        return Y

    def transfer_matrix(self, omega: float, length_scale: float) -> np.ndarray:
        return np.linalg.inv(self.admittance(omega, length_scale))


def tone_metrics(op: LengthOperator, tones, length_scale: float) -> dict:
    z = np.asarray(
        [
            op.transfer_matrix(float(w), length_scale)[
                op.machine.ais, op.machine.input_site
            ]
            for w in tones
        ],
        dtype=complex,
    )
    power = np.abs(z) ** 2
    target = float(power[2])
    off = float(power.sum() - power[2])
    total = float(power.sum())
    return {
        "transfer": z,
        "target_power": target,
        "off_target_power": off,
        "total_power": total,
        "target_fraction": float(target / max(total, 1e-300)),
        "target_gain": float(abs(z[2])),
    }


def grounded_score(metrics: dict, noise_floor: float) -> float:
    return float(
        metrics["target_power"]
        / (metrics["off_target_power"] + float(noise_floor))
    )


def local_length_iteration(
    objective,
    start: float = 1.0,
    log_step: float = 0.22,
    min_log_step: float = 1e-4,
    max_iter: int = 120,
    min_scale: float = 0.05,
    max_scale: float = 50.0,
):
    """Three-point local search; fixed point means no better counted neighbour."""
    u = float(np.log(start))
    step = float(log_step)
    lo = float(np.log(min_scale))
    hi = float(np.log(max_scale))
    history = []
    for k in range(int(max_iter)):
        candidates_u = np.asarray(
            [max(lo, u - step), u, min(hi, u + step)], dtype=float
        )
        candidates_s = np.exp(candidates_u)
        values = np.asarray([objective(float(s)) for s in candidates_s], dtype=float)
        best = int(np.argmax(values))
        history.append(
            {
                "iteration": k,
                "length_scale": float(np.exp(u)),
                "step_log": step,
                "candidate_scales": [float(x) for x in candidates_s],
                "candidate_scores": [float(x) for x in values],
                "choice": best - 1,
            }
        )
        if best == 1 or candidates_u[best] == u:
            step *= 0.5
            if step < min_log_step:
                break
        else:
            u = float(candidates_u[best])
    return float(np.exp(u)), history


def run() -> dict:
    machine = DendriteAIS(branch_len=12)
    parent = int(machine.branch[-2])
    site = int(machine.branch[-1])
    op = LengthOperator(machine, parent, site)

    # Choose the target once from the baseline operator, then freeze it.
    omegas = log_frequency_grid(low=0.002, high=0.8, count=420)
    baseline_gains = np.asarray(
        [
            abs(machine.transfer_matrix(w)[machine.ais, machine.input_site])
            for w in omegas
        ],
        dtype=float,
    )
    omega_star = float(omegas[int(np.argmax(baseline_gains))])
    tones = np.asarray(
        [0.0, 0.35 * omega_star, omega_star, 2.5 * omega_star, 6.0 * omega_star],
        dtype=float,
    )
    baseline = tone_metrics(op, tones, 1.0)

    # ---- ATTACK 1: purity alone -------------------------------------------------
    purity_scales = np.geomspace(0.35, 50.0, 260)
    purity_values = np.asarray(
        [tone_metrics(op, tones, s)["target_fraction"] for s in purity_scales]
    )
    purity_best_i = int(np.argmax(purity_values))
    purity_best_scale = float(purity_scales[purity_best_i])
    purity_best_metrics = tone_metrics(op, tones, purity_best_scale)

    purity_objective = lambda s: tone_metrics(op, tones, s)["target_fraction"]
    purity_iter_scale, purity_history = local_length_iteration(
        purity_objective, start=1.0, max_scale=50.0
    )
    purity_iter_metrics = tone_metrics(op, tones, purity_iter_scale)

    # ---- GROUNDED observer ------------------------------------------------------
    # A bounded readout has an absolute floor.  Express it as a fraction of the
    # baseline target power and report several floors rather than hiding the
    # dependence inside one chosen constant.
    noise_fractions = [0.001, 0.003, 0.01, 0.03, 0.1]
    grounded_rows = []
    grounded_scales = np.geomspace(0.20, 20.0, 320)

    for noise_fraction in noise_fractions:
        noise_floor = float(noise_fraction * baseline["target_power"])
        scores = []
        for s in grounded_scales:
            scores.append(grounded_score(tone_metrics(op, tones, s), noise_floor))
        scores = np.asarray(scores, dtype=float)
        best_i = int(np.argmax(scores))
        broad_scale = float(grounded_scales[best_i])

        objective = lambda s, nf=noise_floor: grounded_score(
            tone_metrics(op, tones, s), nf
        )
        fixed_scale, history = local_length_iteration(
            objective,
            start=1.0,
            min_scale=0.20,
            max_scale=20.0,
        )
        m = tone_metrics(op, tones, fixed_scale)
        h = 1e-3
        derivative = float(
            (
                objective(fixed_scale * np.exp(h))
                - objective(fixed_scale * np.exp(-h))
            )
            / (2.0 * h)
        )
        grounded_rows.append(
            {
                "observer_noise_fraction_of_baseline_target_power": noise_fraction,
                "noise_floor_power": noise_floor,
                "broad_best_scale": broad_scale,
                "broad_best_is_interior": bool(0 < best_i < len(scores) - 1),
                "fixed_point_scale": fixed_scale,
                "fixed_point_score": objective(fixed_scale),
                "dscore_dlogL": derivative,
                "target_fraction": m["target_fraction"],
                "target_gain": m["target_gain"],
                "target_gain_vs_baseline": float(
                    m["target_gain"] / baseline["target_gain"]
                ),
                "iterations": len(history),
                "last_steps": history[-5:],
            }
        )

    # Use the 1% floor as a transparent reference case for operator diagnostics.
    reference = grounded_rows[2]
    natural_scale = float(reference["fixed_point_scale"])
    natural = tone_metrics(op, tones, natural_scale)

    H0 = op.transfer_matrix(omega_star, 1.0)
    Hn = op.transfer_matrix(omega_star, natural_scale)
    z0 = complex(H0[machine.ais, machine.input_site])
    zn = complex(Hn[machine.ais, machine.input_site])
    phase_shift = float(np.angle(zn / z0))

    v0 = baseline["transfer"]
    vn = natural["transfer"]
    scalar = np.vdot(v0, vn) / max(float(np.vdot(v0, v0).real), 1e-15)
    scalar_residual = float(np.linalg.norm(vn - scalar * v0) / np.linalg.norm(vn))

    return {
        "model": {
            "parent_site": parent,
            "length_site": site,
            "baseline_axial_conductance": op.g_ax0,
            "baseline_capacitance": op.c0,
            "baseline_leak": op.g_leak0,
            "geometry_rule": (
                "L->sL at fixed radius/density: g_ax->g_ax/s, "
                "C/leak/quasi-active membrane admittance -> s times local value"
            ),
        },
        "fixed_target": {
            "omega_star": omega_star,
            "tones": [float(x) for x in tones],
            "baseline_target_fraction": baseline["target_fraction"],
            "baseline_target_gain": baseline["target_gain"],
            "baseline_target_power": baseline["target_power"],
            "baseline_off_target_power": baseline["off_target_power"],
        },
        "purity_only_attack": {
            "broad_sweep_max_scale": float(purity_scales[-1]),
            "broad_best_scale": purity_best_scale,
            "broad_best_is_at_upper_boundary": bool(
                purity_best_i == len(purity_scales) - 1
            ),
            "best_target_fraction": purity_best_metrics["target_fraction"],
            "best_target_gain_vs_baseline": float(
                purity_best_metrics["target_gain"] / baseline["target_gain"]
            ),
            "local_iteration_scale": purity_iter_scale,
            "local_iteration_hit_upper_bound": bool(purity_iter_scale > 49.9),
            "local_iteration_target_fraction": purity_iter_metrics["target_fraction"],
            "local_iteration_target_gain_vs_baseline": float(
                purity_iter_metrics["target_gain"] / baseline["target_gain"]
            ),
            "interpretation": (
                "Purity alone rewards isolation: the target fraction rises while "
                "absolute transmission collapses. There is no useful finite natural "
                "length under this ungrounded objective."
            ),
            "last_steps": purity_history[-5:],
        },
        "grounded_bounded_observer": {
            "score": "target_power / (off_target_power + observer_noise_floor)",
            "rows": grounded_rows,
            "reference_noise_fraction": 0.01,
            "reference_fixed_point_scale": natural_scale,
            "reference_target_fraction": natural["target_fraction"],
            "reference_target_gain_vs_baseline": float(
                natural["target_gain"] / baseline["target_gain"]
            ),
        },
        "operator_change_at_reference_length": {
            "target_transfer_phase_shift_radians": phase_shift,
            "five_tone_best_scalar_residual": scalar_residual,
            "interpretation": (
                "The finite length edit changes complex frequency response and is not "
                "equivalent to one scalar multiplier across the five tones."
            ),
        },
        "sigh_bridge": (
            "SighImageSuper taught that a surviving eigenmode can be visually dominant "
            "even while absolute amplitude approaches zero. Gate 12 finds the geometric "
            "version: purity alone drives pathological isolation. A bounded observer's "
            "absolute noise floor grounds the length update and can create a finite fixed point."
        ),
        "claim_boundary": (
            "A finite fixed point here is task- and observer-relative, not a universal "
            "biological preferred length. This does not show that dendrites measure phase "
            "error, execute the search rule, or grow to literal cavity wavelengths."
        ),
    }


if __name__ == "__main__":
    result = run()
    print(json.dumps(result, indent=2, allow_nan=False))
    out = Path("results")
    out.mkdir(exist_ok=True)
    (out / "gate12_natural_length.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
