"""Gate 12: a geometric length operator for dendritic mode purification.

The question is intentionally narrower than "do dendrites grow toward an
eigenmode?".  We ask whether a local segment-length parameter creates a
non-trivial fixed point for the *operator* itself: a length at which nearby
length edits no longer improve the fraction of a chosen temporal mode that
reaches the AIS.

This is a compartmental geometry proxy, not a travelling-wave model.  For one
segment we scale length L -> s L while keeping radius and channel densities
fixed.  In the discrete cable this changes several physical terms together:

    axial conductance      g_ax -> g_ax / s
    membrane capacitance  C    -> s C
    leak conductance      g_L  -> s g_L
    quasi-active density  g_q  -> s g_q

Thus a length edit is not a scalar weight edit.  It changes attenuation and the
complex frequency response through the same morphology/membrane operator.

The SighImageSuper analogy is a fixed-point analogy only:

    state iteration:       x* = F(x*)
    length iteration:      L* = G(L*)

where G takes a local length to whichever nearby length gives better target
purity.  A "natural length" here means a local optimum under that explicitly
counted update rule, not a claim that biological dendrites literally execute
this optimizer.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from active_resonance import log_frequency_grid
from gate11_ais_receipt import DendriteAIS, _edge_outer, _target_fraction


class LengthOperator:
    """One-segment length edit applied to the Gate-11 dendrite+AIS operator."""

    def __init__(self, machine: DendriteAIS, parent: int, site: int):
        self.machine = machine
        self.parent = int(parent)
        self.site = int(site)
        if self.site >= machine.n_d or self.parent >= machine.n_d:
            raise ValueError("Gate 12 currently edits a dendritic segment only")

        # In the passive G matrix an axial edge contributes -g off diagonal.
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

        # Replace original edge g_ax0 by g_ax0/s.
        delta_g_ax = self.g_ax0 * (1.0 / s - 1.0)
        Y += delta_g_ax * _edge_outer(
            self.machine.n, self.parent, self.site
        )

        # Membrane area of this compartment scales with length at fixed radius.
        # The base admittance already contains each term once, so add only the
        # difference from s=1.
        q0 = complex(self.machine.dendrite.quasi_admittance(w)[self.site])
        membrane0 = self.g_leak0 + 1j * w * self.c0 + q0
        Y[self.site, self.site] += (s - 1.0) * membrane0
        return Y

    def transfer_matrix(self, omega: float, length_scale: float) -> np.ndarray:
        return np.linalg.inv(self.admittance(omega, length_scale))


def target_fraction(op: LengthOperator, tones, length_scale: float) -> float:
    Hs = [op.transfer_matrix(float(w), length_scale) for w in tones]
    return _target_fraction(Hs, op.machine.ais, op.machine.input_site)


def local_length_iteration(
    op: LengthOperator,
    tones,
    start: float = 1.0,
    log_step: float = 0.22,
    min_log_step: float = 1e-4,
    max_iter: int = 80,
):
    """Three-point local search whose fixed point is "no better nearby length"."""
    u = float(np.log(start))
    step = float(log_step)
    history = []
    for k in range(int(max_iter)):
        candidates_u = np.asarray([u - step, u, u + step], dtype=float)
        candidates_s = np.exp(candidates_u)
        values = np.asarray(
            [target_fraction(op, tones, s) for s in candidates_s], dtype=float
        )
        best = int(np.argmax(values))
        history.append(
            {
                "iteration": k,
                "length_scale": float(np.exp(u)),
                "step_log": step,
                "candidate_scales": [float(x) for x in candidates_s],
                "candidate_purities": [float(x) for x in values],
                "choice": best - 1,
            }
        )
        if best == 1:
            step *= 0.5
            if step < min_log_step:
                break
        else:
            u = float(candidates_u[best])
    return float(np.exp(u)), history


def run() -> dict:
    machine = DendriteAIS(branch_len=12)

    # Use the same branch edge selected by Gate 11's forward-eligibility/receipt
    # overlap.  That gate selected the terminal branch-A edge 11 -> 12.
    parent = int(machine.branch[-2])
    site = int(machine.branch[-1])
    op = LengthOperator(machine, parent, site)

    # Baseline target carrier is chosen from the unedited operator, then held
    # fixed while length changes.  Otherwise the system could "improve" merely
    # by redefining the target after every edit.
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

    baseline = target_fraction(op, tones, 1.0)

    # Broad audit first: is there an interior geometry optimum at all?
    scales = np.geomspace(0.35, 3.0, 241)
    purities = np.asarray([target_fraction(op, tones, s) for s in scales])
    best_i = int(np.argmax(purities))
    best_scale = float(scales[best_i])
    best_purity = float(purities[best_i])

    # Then apply the explicit Sigh-style fixed-point rule: repeatedly ask whether
    # a slightly shorter/current/slightly longer segment is better.
    natural_scale, history = local_length_iteration(op, tones, start=1.0)
    natural_purity = target_fraction(op, tones, natural_scale)

    # Measure stationarity in log-length coordinates and how much the complex
    # target transfer phase moved.  This is cable phase response, not a claim of
    # a literal wavelength-matching standing wave.
    h = 1e-3
    p_minus = target_fraction(op, tones, natural_scale * np.exp(-h))
    p_plus = target_fraction(op, tones, natural_scale * np.exp(h))
    dP_dlogL = float((p_plus - p_minus) / (2.0 * h))

    H0 = op.transfer_matrix(omega_star, 1.0)
    Hn = op.transfer_matrix(omega_star, natural_scale)
    z0 = complex(H0[machine.ais, machine.input_site])
    zn = complex(Hn[machine.ais, machine.input_site])
    phase_shift = float(np.angle(zn / z0))

    # A length edit should not be equivalent to one scalar multiplier across the
    # five frequencies.  Fit the best complex scalar and report residual.
    v0 = np.asarray(
        [
            op.transfer_matrix(w, 1.0)[machine.ais, machine.input_site]
            for w in tones
        ],
        dtype=complex,
    )
    vn = np.asarray(
        [
            op.transfer_matrix(w, natural_scale)[machine.ais, machine.input_site]
            for w in tones
        ],
        dtype=complex,
    )
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
            "baseline_target_fraction": baseline,
        },
        "broad_length_sweep": {
            "minimum_scale": float(scales[0]),
            "maximum_scale": float(scales[-1]),
            "best_scale": best_scale,
            "best_target_fraction": best_purity,
            "improvement_over_baseline": float(best_purity - baseline),
            "best_is_interior": bool(0 < best_i < len(scales) - 1),
        },
        "natural_length_iteration": {
            "start_scale": 1.0,
            "final_scale": natural_scale,
            "final_target_fraction": natural_purity,
            "improvement_over_baseline": float(natural_purity - baseline),
            "dP_dlogL_at_final": dP_dlogL,
            "iterations": len(history),
            "last_steps": history[-8:],
        },
        "operator_change": {
            "target_transfer_phase_shift_radians": phase_shift,
            "five_tone_best_scalar_residual": scalar_residual,
            "interpretation": (
                "A segment-length edit changes the frequency-dependent operator; "
                "the five-tone response is not generally reproducible by one scalar weight."
            ),
        },
        "claim_boundary": (
            "This gate asks whether a local compartment-length parameter has an "
            "operator-level fixed point under an explicit local search rule. It is "
            "not evidence that dendrites measure phase error, execute this optimizer, "
            "or grow to literal acoustic wavelengths."
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
