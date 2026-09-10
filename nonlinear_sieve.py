"""Gate 6 machinery: let the resonant branch select the mode *before* nonlinearity.

This module adds one bounded regenerative hotspot to the quasi-active dendrite.
The nonlinearity is intentionally generic. It is not claimed to be an NMDA
kinetic model. The question is architectural:

    distributed resonant operator -> local mode concentration -> nonlinear event

If the branch is doing useful mode purification, equal-budget inputs that are
phase/frequency matched to the branch should cross the local nonlinear regime
more strongly than phase-scrambled, off-band, or spatially separated controls.
"""

from __future__ import annotations

import numpy as np

from active_resonance import QuasiActiveDendrite, log_frequency_grid, peak_info


class ResonantNonlinearSieve:
    def __init__(self, model: QuasiActiveDendrite | None = None):
        self.model = model or QuasiActiveDendrite()
        self.A = self.model.state_matrix()
        self.nv = self.model.n
        self.ns = self.A.shape[0]

        # Keep the nonlinear site proximal enough to be a meeting zone while
        # still inside branch A rather than at the soma.
        self.hotspot = int(self.model.branch_a[2])
        self.site_a = int(self.model.branch_a[-1])
        self.site_b = int(self.model.branch_a[-5])
        self.separated_site_b = int(self.model.branch_b[-1])

        ws = log_frequency_grid(0.002, 1.0, 600)
        g = self.model.spectrum(ws, self.site_a, self.model.soma)
        self.omega = float(peak_info(ws, g)["omega"])

        # Calibrate the two same-branch sources so their *linear* phasors have
        # equal amplitude at the hotspot at the branch's preferred frequency.
        H = self.model.transfer_matrix(self.omega)
        ha = complex(H[self.hotspot, self.site_a])
        hb = complex(H[self.hotspot, self.site_b])
        self.rel_amp_a = 1.0
        self.rel_amp_b = float(abs(ha) / abs(hb))
        self.phase_a = 0.0
        self.phase_b_matched = float(np.angle(ha) - np.angle(hb))

        # Set the matched *linear* hotspot oscillation to a fixed amplitude.
        # This is an operating-point calibration, not learned fitting.
        desired_hotspot_amp = 0.080
        linear_sum = abs(ha) * self.rel_amp_a + abs(hb) * self.rel_amp_b
        self.drive_scale = float(desired_hotspot_amp / linear_sum)

        # Singles are ~0.04, quadrature pair ~0.0566, matched pair ~0.08 in
        # the linear phasor prediction. Put the regenerative knee between
        # quadrature and matched. The current is bounded, so the ODE remains
        # dissipative at large voltage even though the local feedback is positive.
        self.threshold = 0.060
        self.slope = 0.003
        self.i_max = 0.018
        self.v_sat = 0.050

    def same_protocol_on(self, model: QuasiActiveDendrite):
        """Return another substrate driven by the *identical* Gate-6 protocol.

        This is used for the Gate-6b attacker. Nothing is recalibrated after the
        substrate is changed: carrier, phases, source amplitudes, nonlinear knee,
        and source/hotspot addresses are copied from the active machine.
        """
        other = object.__new__(ResonantNonlinearSieve)
        other.model = model
        other.A = model.state_matrix()
        other.nv = model.n
        other.ns = other.A.shape[0]
        for name in (
            "hotspot", "site_a", "site_b", "separated_site_b", "omega",
            "rel_amp_a", "rel_amp_b", "phase_a", "phase_b_matched",
            "drive_scale", "threshold", "slope", "i_max", "v_sat",
        ):
            setattr(other, name, getattr(self, name))
        return other

    def input_vector(self, site: int) -> np.ndarray:
        b = np.zeros(self.ns, dtype=float)
        b[int(site)] = 1.0 / self.model.C[int(site)]
        return b

    def regenerative_current(self, v: float) -> float:
        a = abs(float(v))
        z = np.clip((a - self.threshold) / self.slope, -60.0, 60.0)
        gate = 1.0 / (1.0 + np.exp(-z))
        return float(self.i_max * gate * np.tanh(float(v) / self.v_sat))

    def phase_to_align(self, omega: float, site: int) -> float:
        H = self.model.transfer_matrix(float(omega))
        ha = complex(H[self.hotspot, self.site_a])
        hx = complex(H[self.hotspot, int(site)])
        return float(np.angle(ha) - np.angle(hx))

    def case_specs(self) -> dict:
        w = self.omega
        low = 0.35 * w
        high = 2.5 * w
        return {
            "matched": {
                "site_b": self.site_b,
                "omega": w,
                "phase_b": self.phase_b_matched,
                "use_a": True,
                "use_b": True,
            },
            "quadrature": {
                "site_b": self.site_b,
                "omega": w,
                "phase_b": self.phase_b_matched + np.pi / 2.0,
                "use_a": True,
                "use_b": True,
            },
            "antiphase": {
                "site_b": self.site_b,
                "omega": w,
                "phase_b": self.phase_b_matched + np.pi,
                "use_a": True,
                "use_b": True,
            },
            "wrong_low_carrier": {
                "site_b": self.site_b,
                "omega": low,
                "phase_b": self.phase_to_align(low, self.site_b),
                "use_a": True,
                "use_b": True,
            },
            "wrong_high_carrier": {
                "site_b": self.site_b,
                "omega": high,
                "phase_b": self.phase_to_align(high, self.site_b),
                "use_a": True,
                "use_b": True,
            },
            "spatially_separated": {
                "site_b": self.separated_site_b,
                "omega": w,
                "phase_b": self.phase_to_align(w, self.separated_site_b),
                "use_a": True,
                "use_b": True,
            },
            "single_A": {
                "site_b": self.site_b,
                "omega": w,
                "phase_b": self.phase_b_matched,
                "use_a": True,
                "use_b": False,
            },
            "single_B": {
                "site_b": self.site_b,
                "omega": w,
                "phase_b": self.phase_b_matched,
                "use_a": False,
                "use_b": True,
            },
        }

    def _rhs(self, t: float, x: np.ndarray, spec: dict, nonlinear: bool) -> np.ndarray:
        y = self.A @ x
        w = float(spec["omega"])
        if spec["use_a"]:
            ia = self.drive_scale * self.rel_amp_a * np.sin(w * t + self.phase_a)
            y += ia * self.input_vector(self.site_a)
        if spec["use_b"]:
            ib = self.drive_scale * self.rel_amp_b * np.sin(w * t + float(spec["phase_b"]))
            y += ib * self.input_vector(int(spec["site_b"]))
        if nonlinear:
            ireg = self.regenerative_current(x[self.hotspot])
            y[self.hotspot] += ireg / self.model.C[self.hotspot]
        return y

    def simulate(self, spec: dict, nonlinear: bool, periods: int = 22, keep_periods: int = 8, dt: float = 0.2):
        w = float(spec["omega"])
        period = 2.0 * np.pi / w
        total_steps = int(np.ceil(periods * period / dt))
        keep_steps = int(np.ceil(keep_periods * period / dt))
        start_keep = max(0, total_steps - keep_steps)

        x = np.zeros(self.ns, dtype=float)
        soma = []
        hot = []
        for step in range(total_steps):
            t = step * dt
            k1 = self._rhs(t, x, spec, nonlinear)
            k2 = self._rhs(t + 0.5 * dt, x + 0.5 * dt * k1, spec, nonlinear)
            k3 = self._rhs(t + 0.5 * dt, x + 0.5 * dt * k2, spec, nonlinear)
            k4 = self._rhs(t + dt, x + dt * k3, spec, nonlinear)
            x = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            if step >= start_keep:
                soma.append(float(x[self.model.soma]))
                hot.append(float(x[self.hotspot]))
        return np.asarray(soma), np.asarray(hot)

    def evaluate_case(self, spec: dict) -> dict:
        soma_lin, hot_lin = self.simulate(spec, nonlinear=False)
        soma_nl, hot_nl = self.simulate(spec, nonlinear=True)
        residual = soma_nl - soma_lin

        def rms(x):
            return float(np.sqrt(np.mean(np.square(x))))

        return {
            "omega": float(spec["omega"]),
            "site_b": int(spec["site_b"]),
            "phase_b": float(spec["phase_b"]),
            "linear_hotspot_peak": float(np.max(np.abs(hot_lin))),
            "linear_hotspot_threshold_occupancy": float(np.mean(np.abs(hot_lin) > self.threshold)),
            "nonlinear_hotspot_peak": float(np.max(np.abs(hot_nl))),
            "linear_soma_rms": rms(soma_lin),
            "nonlinear_soma_rms": rms(soma_nl),
            "nonlinear_residual_soma_rms": rms(residual),
            "residual_over_linear_soma": rms(residual) / max(rms(soma_lin), 1e-15),
        }

    def run_all(self) -> dict:
        specs = self.case_specs()
        return {name: self.evaluate_case(spec) for name, spec in specs.items()}
