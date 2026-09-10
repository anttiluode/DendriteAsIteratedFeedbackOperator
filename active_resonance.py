"""Quasi-active dendritic resonance as an explicit mode-selective operator.

This is deliberately a *phenomenological* restorative membrane model, not a
claim that a biological dendrite is a lossless waveguide.  Each active
compartment gets one slow gating variable z.  In the linearized frequency
domain the gate contributes

    g_q / (1 + i omega tau_q)

so slow/DC voltage changes are strongly opposed, while at intermediate
frequency the restorative gate cannot fully follow and the membrane can show
a genuine band-pass impedance peak.

The important object is still the operator.  At each temporal frequency,

    H(omega) = [G + i omega C + D_q(omega)]^-1,

and local edits to D_q change the global transfer function.
"""

from __future__ import annotations

import numpy as np

from dendrite_operator import BranchedCable


class QuasiActiveDendrite:
    """Branched cable with one restorative quasi-active gate per dendritic compartment.

    The two branches intentionally have different gate time constants so that
    the same morphology supports two distinct non-zero temporal pass bands.
    This is an operator laboratory, not a fitted ion-channel model.
    """

    def __init__(
        self,
        branch_len: int = 12,
        g_a: float = 0.30,
        tau_a: float = 18.0,
        g_b: float = 0.18,
        tau_b: float = 42.0,
    ):
        self.base = BranchedCable(branch_len=branch_len)
        self.n = self.base.n
        self.soma = self.base.soma
        self.branch_a = self.base.branch_a
        self.branch_b = self.base.branch_b
        self.C = self.base.C.copy()
        self.G = self.base.G.copy()

        self.g_q = np.zeros(self.n, dtype=float)
        self.tau_q = np.ones(self.n, dtype=float)
        self.g_q[self.branch_a] = float(g_a)
        self.g_q[self.branch_b] = float(g_b)
        self.tau_q[self.branch_a] = float(tau_a)
        self.tau_q[self.branch_b] = float(tau_b)
        self.active_sites = np.concatenate([self.branch_a, self.branch_b]).astype(int)

    def quasi_admittance(self, omega: float) -> np.ndarray:
        """Return the complex quasi-active membrane admittance per compartment."""
        w = float(omega)
        return self.g_q / (1.0 + 1j * w * self.tau_q)

    def admittance_matrix(self, omega: float) -> np.ndarray:
        """Voltage-domain operator Y(omega) whose inverse is the transfer matrix."""
        w = float(omega)
        return (
            self.G.astype(complex)
            + 1j * w * np.diag(self.C)
            + np.diag(self.quasi_admittance(w))
        )

    def transfer_matrix(self, omega: float) -> np.ndarray:
        """Current-to-voltage transfer H(omega)=Y(omega)^-1."""
        return np.linalg.inv(self.admittance_matrix(omega))

    def transfer(self, omega: float, input_site: int, output_site: int | None = None) -> complex:
        if output_site is None:
            output_site = self.soma
        H = self.transfer_matrix(omega)
        return complex(H[int(output_site), int(input_site)])

    def spectrum(self, omegas, input_site: int, output_site: int | None = None) -> np.ndarray:
        return np.asarray(
            [abs(self.transfer(w, input_site, output_site)) for w in omegas],
            dtype=float,
        )

    def passive_spectrum(self, omegas, input_site: int, output_site: int | None = None) -> np.ndarray:
        """Matched morphology with the quasi-active gates removed."""
        if output_site is None:
            output_site = self.soma
        out = []
        Cdiag = np.diag(self.C)
        for w in omegas:
            H = np.linalg.inv(self.G.astype(complex) + 1j * float(w) * Cdiag)
            out.append(abs(H[int(output_site), int(input_site)]))
        return np.asarray(out, dtype=float)

    def state_matrix(self) -> np.ndarray:
        """Continuous-time real state matrix for [V, z_active].

        C dV/dt = -G V - R z + I
        tau dz/dt = V_site - z

        Eliminating z under sinusoidal drive gives the frequency-domain
        quasi-admittance used above.
        """
        sites = self.active_sites
        m = len(sites)
        cinv = 1.0 / self.C

        A_vv = -(cinv[:, None] * self.G)
        A_vz = np.zeros((self.n, m), dtype=float)
        A_zv = np.zeros((m, self.n), dtype=float)
        A_zz = np.zeros((m, m), dtype=float)

        for k, site in enumerate(sites):
            tau = self.tau_q[site]
            A_vz[site, k] = -self.g_q[site] / self.C[site]
            A_zv[k, site] = 1.0 / tau
            A_zz[k, k] = -1.0 / tau

        return np.block([[A_vv, A_vz], [A_zv, A_zz]])

    def poles(self) -> np.ndarray:
        return np.linalg.eigvals(self.state_matrix())

    def edited_transfer_matrix(self, omega: float, site: int, delta_g: float) -> np.ndarray:
        """Transfer after one local quasi-active conductance edit.

        At fixed omega this is a rank-one diagonal admittance edit, so the
        corresponding global resolvent change should be rank one up to roundoff.
        """
        site = int(site)
        if self.g_q[site] == 0.0:
            raise ValueError("edited site must carry a quasi-active gate")
        Y = self.admittance_matrix(omega).copy()
        Y[site, site] += float(delta_g) / (1.0 + 1j * float(omega) * self.tau_q[site])
        return np.linalg.inv(Y)


def log_frequency_grid(low: float = 0.002, high: float = 1.0, count: int = 500) -> np.ndarray:
    """Frequency grid including DC followed by logarithmically spaced omega > 0."""
    return np.concatenate([[0.0], np.geomspace(float(low), float(high), int(count))])


def peak_info(omegas, gains) -> dict:
    omegas = np.asarray(omegas, dtype=float)
    gains = np.asarray(gains, dtype=float)
    idx = int(np.argmax(gains))
    peak = float(gains[idx])
    half_power_gain = peak / np.sqrt(2.0)

    left = None
    for i in range(idx - 1, -1, -1):
        if gains[i] <= half_power_gain:
            left = float(omegas[i])
            break
    right = None
    for i in range(idx + 1, len(gains)):
        if gains[i] <= half_power_gain:
            right = float(omegas[i])
            break

    bandwidth = None
    q = None
    if left is not None and right is not None and right > left and omegas[idx] > 0:
        bandwidth = right - left
        q = float(omegas[idx] / bandwidth)

    return {
        "index": idx,
        "omega": float(omegas[idx]),
        "gain": peak,
        "dc_gain": float(gains[0]),
        "peak_over_dc": float(peak / gains[0]) if gains[0] > 0 else float("inf"),
        "half_power_left": left,
        "half_power_right": right,
        "bandwidth": bandwidth,
        "q_like": q,
    }


def equal_tone_purification(model: QuasiActiveDendrite, input_site: int, omega_target: float, output_site=None):
    """Operational 'purifier' metric for five equal-amplitude input tones.

    Input energy is 20% in each tone.  The output target fraction is computed
    only from linear transfer gains; no information is created.  A value above
    0.2 means the dendritic operator enriched the selected temporal mode.
    """
    w = float(omega_target)
    tones = np.asarray([0.0, 0.35 * w, w, 2.5 * w, 6.0 * w], dtype=float)
    active = model.spectrum(tones, input_site, output_site)
    passive = model.passive_spectrum(tones, input_site, output_site)

    def frac(g):
        p = np.square(g)
        return float(p[2] / p.sum())

    fa = frac(active)
    fp = frac(passive)
    return {
        "tones": [float(x) for x in tones],
        "input_target_energy_fraction": 0.2,
        "active_output_target_energy_fraction": fa,
        "passive_output_target_energy_fraction": fp,
        "active_enrichment_over_input": fa / 0.2,
        "active_enrichment_over_passive": fa / fp if fp > 0 else float("inf"),
        "active_gains": [float(x) for x in active],
        "passive_gains": [float(x) for x in passive],
    }
