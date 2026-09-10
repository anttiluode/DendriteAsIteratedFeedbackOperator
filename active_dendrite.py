import numpy as np


def resonant_state_matrix(cable, sites, g_r=0.2, tau=20.0):
    """Linear delayed-restoring current on selected dendritic compartments.

    This is a deliberately small subthreshold resonator:
        C dV/dt = -G V - g_r S w + I
        tau dw/dt = S^T V - w

    It is not claimed as a quantitative HCN model.  It is the minimal
    established dynamical motif: a delayed restorative current can turn an
    otherwise low-pass dendritic cable into a resonant transfer system.
    """
    sites = np.asarray(sites, dtype=int)
    n = cable.n
    m = len(sites)
    A = np.zeros((n + m, n + m), dtype=float)
    A[:n, :n] = -(cable.G / cable.C[:, None])

    for k, site in enumerate(sites):
        A[site, n + k] = -float(g_r) / cable.C[site]
        A[n + k, site] = 1.0 / float(tau)
        A[n + k, n + k] = -1.0 / float(tau)
    return A


def frequency_response(cable, A, input_site, output_site, omegas):
    """Complex transfer from a unit current at input_site to V(output_site)."""
    omegas = np.asarray(omegas, dtype=float)
    B = np.zeros(A.shape[0], dtype=float)
    B[int(input_site)] = 1.0 / cable.C[int(input_site)]
    C = np.zeros(A.shape[0], dtype=float)
    C[int(output_site)] = 1.0
    eye = np.eye(A.shape[0])
    out = np.empty(len(omegas), dtype=complex)
    for k, omega in enumerate(omegas):
        out[k] = C @ np.linalg.solve(1j * omega * eye - A, B)
    return out


def nmda_block(v, vhalf=0.30, slope=0.08):
    """Smooth voltage-dependent gate used for a qualitative NMDA-like toy.

    Voltage is dimensionless in this laboratory.  The function is a sigmoid
    approximation to voltage-dependent relief of block; it is not a fitted
    receptor model.
    """
    z = np.clip(-(np.asarray(v) - float(vhalf)) / float(slope), -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(z))


def nmda_current_slope(v, vhalf=0.30, slope=0.08, reversal=1.0):
    """d/dV [ B(V) (E - V) ]."""
    b = nmda_block(v, vhalf=vhalf, slope=slope)
    bp = b * (1.0 - b) / float(slope)
    return bp * (float(reversal) - v) - b


def nmda_rhs(
    cable,
    v,
    t,
    events,
    g_peak=0.4,
    tau_s=8.0,
    vhalf=0.30,
    slope=0.08,
    reversal=1.0,
):
    """Voltage derivative for passive cable plus local NMDA-like conductances."""
    v = np.asarray(v, dtype=float)
    current = -cable.G @ v
    for site, t0 in events:
        if t >= t0:
            syn = np.exp(-(t - t0) / float(tau_s))
            block = nmda_block(v[int(site)], vhalf=vhalf, slope=slope)
            current[int(site)] += (
                float(g_peak)
                * syn
                * block
                * (float(reversal) - v[int(site)])
            )
    return current / cable.C


def nmda_jacobian(
    cable,
    v,
    t,
    events,
    g_peak=0.4,
    tau_s=8.0,
    vhalf=0.30,
    slope=0.08,
    reversal=1.0,
):
    """Instantaneous Jacobian d(dV/dt)/dV.

    The external synaptic gates are held fixed.  Any difference between two
    Jacobians evaluated with the same events/time but different V therefore
    comes from voltage-dependent state, not from a different input schedule.
    """
    J = -(cable.G / cable.C[:, None]).copy()
    for site, t0 in events:
        site = int(site)
        if t >= t0:
            syn = np.exp(-(t - t0) / float(tau_s))
            slope_i = nmda_current_slope(
                v[site], vhalf=vhalf, slope=slope, reversal=reversal
            )
            J[site, site] += float(g_peak) * syn * slope_i / cable.C[site]
    return J


def simulate_nmda(
    cable,
    events,
    g_peak=0.4,
    tau_s=8.0,
    vhalf=0.30,
    slope=0.08,
    reversal=1.0,
    dt=0.05,
    duration=60.0,
):
    """Deterministic RK4 simulation. Returns (times, voltage_history)."""
    times = np.arange(0.0, float(duration), float(dt))
    history = np.empty((len(times), cable.n), dtype=float)
    v = np.zeros(cable.n, dtype=float)

    def f(vv, tt):
        return nmda_rhs(
            cable,
            vv,
            tt,
            events,
            g_peak=g_peak,
            tau_s=tau_s,
            vhalf=vhalf,
            slope=slope,
            reversal=reversal,
        )

    for k, t in enumerate(times):
        history[k] = v
        h = float(dt)
        k1 = f(v, t)
        k2 = f(v + 0.5 * h * k1, t + 0.5 * h)
        k3 = f(v + 0.5 * h * k2, t + 0.5 * h)
        k4 = f(v + h * k3, t + h)
        v = v + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    return times, history
