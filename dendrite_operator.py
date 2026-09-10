import numpy as np


def _add_edge(G, i, j, g):
    b = np.zeros(G.shape[0])
    b[i] = 1.0
    b[j] = -1.0
    G += g * np.outer(b, b)


class BranchedCable:
    """Small passive compartmental cable used as an operator laboratory.

    Dynamics: C dV/dt = -G V + I(t)
    where G contains leak + axial conductances.  C is diagonal.
    """

    def __init__(self, branch_len=12):
        self.branch_len = int(branch_len)
        self.n = 1 + 2 * self.branch_len
        self.soma = 0
        self.branch_a = np.arange(1, 1 + self.branch_len)
        self.branch_b = np.arange(1 + self.branch_len, self.n)

        # Capacitance is a rough compartment volume/area proxy.  The soma is larger.
        self.C = np.ones(self.n)
        self.C[self.soma] = 5.0
        self.C[self.branch_a[-1]] = 2.2
        self.C[self.branch_b[-1]] = 2.2

        # Leak makes every autonomous mode decay.  A slight branch asymmetry is
        # deliberate: biology is not a perfectly symmetric Y and scalar weights
        # cannot capture the resulting temporal kernels.
        self.g_leak = np.full(self.n, 0.055)
        self.g_leak[self.soma] = 0.09
        self.g_leak[self.branch_a] *= 0.95
        self.g_leak[self.branch_b] *= 1.08

        self.G = np.diag(self.g_leak.copy())

        # Soma-to-branch necks.
        _add_edge(self.G, self.soma, self.branch_a[0], 0.38)
        _add_edge(self.G, self.soma, self.branch_b[0], 0.22)

        # Branch A: relatively uniform cable.
        for k in range(self.branch_len - 1):
            _add_edge(self.G, self.branch_a[k], self.branch_a[k + 1], 0.72)

        # Branch B: a taper/bottleneck near the middle.
        for k in range(self.branch_len - 1):
            g = 0.68 if k < self.branch_len // 2 else 0.34
            _add_edge(self.G, self.branch_b[k], self.branch_b[k + 1], g)

    def symmetric_generator(self):
        """Return B = C^-1/2 G C^-1/2, similar to the physical generator C^-1 G."""
        c_isqrt = 1.0 / np.sqrt(self.C)
        return (c_isqrt[:, None] * self.G) * c_isqrt[None, :]

    def modes(self):
        B = self.symmetric_generator()
        lam, Q = np.linalg.eigh(B)
        return lam, Q

    def step_operator(self, dt=0.5):
        """Exact passive step P = exp(-dt C^-1 G), constructed by similarity transform."""
        lam, Q = self.modes()
        E = (Q * np.exp(-dt * lam)) @ Q.T
        c_sqrt = np.sqrt(self.C)
        c_isqrt = 1.0 / c_sqrt
        # V = C^-1/2 y, y0 = C^1/2 V0
        return (c_isqrt[:, None] * E) * c_sqrt[None, :]

    def resolvent(self, s=0.15):
        """Continuous-time transfer H(s) = (G + s C)^-1."""
        return np.linalg.inv(self.G + s * np.diag(self.C))

    def edited_edge_resolvent(self, i, j, delta_g, s=0.15):
        G2 = self.G.copy()
        _add_edge(G2, int(i), int(j), float(delta_g))
        return np.linalg.inv(G2 + s * np.diag(self.C))

    def impulse_kernel(self, site, steps=100, dt=0.5):
        P = self.step_operator(dt=dt)
        x = np.zeros(self.n)
        # Unit charge: voltage jump is Q/C locally.
        x[int(site)] = 1.0 / self.C[int(site)]
        out = np.empty(steps)
        for t in range(steps):
            out[t] = x[self.soma]
            x = P @ x
        return out

    def grounded_fixed_point(self, cue, alpha=0.08, dt=0.5):
        P = self.step_operator(dt=dt)
        I = np.eye(self.n)
        cue = np.asarray(cue, dtype=float)
        return np.linalg.solve(I - (1.0 - alpha) * P, alpha * cue)


def effective_rank(svals, eps=1e-15):
    s = np.asarray(svals, dtype=float)
    s = s[s > eps]
    if len(s) == 0:
        return 0.0
    p = s / s.sum()
    return float(np.exp(-(p * np.log(p)).sum()))


def cosine(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))
