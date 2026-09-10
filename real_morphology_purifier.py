"""Real-morphology Gate 10: does the eigenmode-purifier picture survive a real arbor?

This gate reuses the exact pinned human L2/3 morphology already audited in
anttiluode/Operaattori (cell 1125 from the FCI repository).  It does *not* reuse
the Y toy geometry.

The morphology is turned into a sparse compartmental admittance graph using
standard cylinder scalings:

    g_axial = pi r^2 / (R_a L)
    C       = c_m * area
    g_leak  = g_m * area

A deliberately simple restorative quasi-active membrane is then added on
neuritic compartments:

    Y_q(w) = g_q / (1 + i w tau_q)

The scientific question is narrow: with one frozen membrane parameter set,
do several *different real dendritic paths* acquire non-zero transfer peaks and
enrich those preferred temporal components relative to the same passive
morphology?

This is still phenomenological.  It is not a fitted HCN distribution and it is
not evidence that real dendrites attain the measured selectivity in vivo.
"""

from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu


SOURCE_COMMIT = "75ad8b4d81a7f51bf888b30650c543592340db06"
SOURCE_NAME = "2013_03_06_cell11_1125_H41_06.asc"
SOURCE_REL = (
    "simulating_neurons/neuron_models/human/eyal/"
    "Human_L23_PC_0603_11_937_Eyal_passive_dends_simple_soma/"
    "morphologies/" + SOURCE_NAME
)
SOURCE_URL = "https://raw.githubusercontent.com/ido4848/FCI/" + SOURCE_COMMIT + "/" + SOURCE_REL


@dataclass
class Tree:
    positions: np.ndarray
    parents: np.ndarray
    radii_um: np.ndarray
    section_types: np.ndarray


def _enum_int(x) -> int:
    try:
        return int(x)
    except (TypeError, ValueError):
        return int(getattr(x, "value", -1))


def download_source(dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 100_000:
        return dest
    req = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "DendritePurifier/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r, dest.open("wb") as f:
        f.write(r.read())
    return dest


def load_tree(path: Path, duplicate_tol: float = 1e-5) -> Tree:
    from morphio import Morphology

    morph = Morphology(str(path))
    soma_pts = np.asarray(morph.soma.points, dtype=float)
    soma_diams = np.asarray(morph.soma.diameters, dtype=float)
    if len(soma_pts):
        root = np.mean(soma_pts, axis=0)
        root_r = float(np.mean(soma_diams) / 2.0) if len(soma_diams) else 5.0
    else:
        first = [np.asarray(s.points[0], dtype=float) for s in morph.root_sections if len(s.points)]
        root = np.mean(np.stack(first), axis=0)
        root_r = 5.0

    pos = [root]
    parents = [-1]
    radii = [root_r]
    stypes = [1]
    section_end: dict[int, int] = {}

    for section in morph.iter():
        sid = int(section.id)
        stype = _enum_int(section.type)
        if section.is_root:
            parent = 0
        else:
            parent = section_end[int(section.parent.id)]

        pts = np.asarray(section.points, dtype=float)
        diams = np.asarray(section.diameters, dtype=float)
        start = 1 if len(pts) and np.linalg.norm(pts[0] - pos[parent]) <= duplicate_tol else 0
        for j in range(start, len(pts)):
            p = np.asarray(pts[j], dtype=float)
            if np.linalg.norm(p - pos[parent]) <= duplicate_tol:
                continue
            idx = len(pos)
            pos.append(p)
            parents.append(parent)
            radii.append(max(float(diams[j]) / 2.0, 0.05))
            stypes.append(stype)
            parent = idx
        section_end[sid] = parent

    return Tree(
        positions=np.asarray(pos, dtype=float),
        parents=np.asarray(parents, dtype=np.int64),
        radii_um=np.asarray(radii, dtype=float),
        section_types=np.asarray(stypes, dtype=np.int64),
    )


class RealPurifier:
    """Sparse nodal frequency-domain model on the real morphology."""

    def __init__(
        self,
        tree: Tree,
        ra_ohm_cm: float = 150.0,
        cm_f_per_cm2: float = 1.0e-6,
        gm_s_per_cm2: float = 5.0e-5,
        q_over_leak: float = 4.0,
        tau_q_s: float = 0.030,
    ):
        self.tree = tree
        self.n = len(tree.parents)
        self.ra = float(ra_ohm_cm)
        self.cm = float(cm_f_per_cm2)
        self.gm = float(gm_s_per_cm2)
        self.q_over_leak = float(q_over_leak)
        self.tau = float(tau_q_s)

        rows, cols, vals = [], [], []
        diag_ax = np.zeros(self.n, dtype=float)
        area_cm2 = np.zeros(self.n, dtype=float)
        root_r_cm = tree.radii_um[0] * 1e-4
        area_cm2[0] = 4.0 * np.pi * root_r_cm * root_r_cm

        self.edge_length_um = np.zeros(self.n, dtype=float)
        for i in range(1, self.n):
            p = int(tree.parents[i])
            L_um = max(float(np.linalg.norm(tree.positions[i] - tree.positions[p])), 0.05)
            r_um = max(float(tree.radii_um[i]), 0.05)
            self.edge_length_um[i] = L_um
            L_cm = L_um * 1e-4
            r_cm = r_um * 1e-4
            gax = np.pi * r_cm * r_cm / (self.ra * L_cm)
            diag_ax[i] += gax
            diag_ax[p] += gax
            rows.extend([i, p])
            cols.extend([p, i])
            vals.extend([-gax, -gax])
            area_cm2[i] += 2.0 * np.pi * r_cm * L_cm

        rows.extend(range(self.n))
        cols.extend(range(self.n))
        vals.extend(diag_ax.tolist())
        self.Gax = sparse.csr_matrix((vals, (rows, cols)), shape=(self.n, self.n))

        # Tiny floor prevents pathological zero-area synthetic/root cases while
        # remaining negligible relative to ordinary neurite membrane area.
        area_floor = max(float(np.median(area_cm2[area_cm2 > 0])) * 1e-4, 1e-16)
        area_cm2 = np.maximum(area_cm2, area_floor)
        self.area_cm2 = area_cm2
        self.C = self.cm * area_cm2
        self.g_leak = self.gm * area_cm2

        # Conventional morphology type 2 is axon; keep it passive.  Types 3/4
        # (basal/apical dendrite) and any non-axon neurite get the same frozen
        # restorative density.  The soma root is passive in this first gate.
        dendritic = (tree.section_types != 2)
        dendritic[0] = False
        self.g_q = np.where(dendritic, self.q_over_leak * self.g_leak, 0.0)

        child_count = np.zeros(self.n, dtype=np.int64)
        for i in range(1, self.n):
            child_count[int(tree.parents[i])] += 1
        self.child_count = child_count
        self.path_distance_um = np.zeros(self.n, dtype=float)
        for i in range(1, self.n):
            self.path_distance_um[i] = self.path_distance_um[int(tree.parents[i])] + self.edge_length_um[i]

    def matrix(self, frequency_hz: float, active: bool) -> sparse.csc_matrix:
        w = 2.0 * np.pi * float(frequency_hz)
        diag = self.g_leak.astype(np.complex128) + 1j * w * self.C
        if active:
            diag = diag + self.g_q / (1.0 + 1j * w * self.tau)
        return (self.Gax.astype(np.complex128) + sparse.diags(diag)).tocsc()

    def solve_inputs(self, frequency_hz: float, inputs: list[int], active: bool) -> np.ndarray:
        lu = splu(self.matrix(frequency_hz, active))
        B = np.zeros((self.n, len(inputs)), dtype=np.complex128)
        for j, node in enumerate(inputs):
            B[int(node), j] = 1.0
        return np.asarray(lu.solve(B), dtype=np.complex128)

    def select_diverse_tips(self, count: int = 6) -> list[int]:
        tips = np.where((self.child_count == 0) & (self.tree.section_types != 2))[0]
        # Prefer long electrotonic routes while retaining both major dendrite
        # classes when available.
        by_type: dict[int, list[int]] = {}
        for i in tips:
            by_type.setdefault(int(self.tree.section_types[i]), []).append(int(i))
        for t in by_type:
            by_type[t].sort(key=lambda i: self.path_distance_um[i], reverse=True)

        chosen: list[int] = []
        for t in (4, 3):
            chosen.extend(by_type.get(t, [])[: max(1, count // 2)])
        remaining = sorted((int(i) for i in tips if int(i) not in chosen), key=lambda i: self.path_distance_um[i], reverse=True)
        chosen.extend(remaining[: max(0, count - len(chosen))])
        return chosen[:count]

    def proximal_node(self, tip: int) -> int:
        i = int(tip)
        prev = i
        while int(self.tree.parents[i]) > 0:
            prev = i
            i = int(self.tree.parents[i])
        # i is first child of root if parent[i] == 0; `prev` only differs when
        # a malformed path terminates earlier.
        return int(i if int(self.tree.parents[i]) == 0 else prev)


def target_fraction(gains: np.ndarray, target_index: int = 2) -> float:
    p = np.square(np.asarray(gains, dtype=float))
    return float(p[target_index] / np.sum(p))


def run(out_dir: Path | str = "results/real_morphology") -> dict:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    source = download_source(out / "_source" / SOURCE_NAME)
    tree = load_tree(source)
    model = RealPurifier(tree)
    tips = model.select_diverse_tips(6)

    freqs = np.concatenate([[0.0], np.geomspace(0.5, 100.0, 72)])
    active_soma = np.zeros((len(freqs), len(tips)), dtype=float)
    passive_soma = np.zeros_like(active_soma)
    for k, f in enumerate(freqs):
        Xa = model.solve_inputs(float(f), tips, active=True)
        Xp = model.solve_inputs(float(f), tips, active=False)
        active_soma[k] = np.abs(Xa[0, :])
        passive_soma[k] = np.abs(Xp[0, :])

    rows = []
    for j, tip in enumerate(tips):
        ia = int(np.argmax(active_soma[:, j]))
        ip = int(np.argmax(passive_soma[:, j]))
        fpeak = float(freqs[ia])
        prox = model.proximal_node(tip)
        if fpeak <= 0.0:
            rows.append({
                "tip": int(tip),
                "section_type": int(tree.section_types[tip]),
                "path_distance_um": float(model.path_distance_um[tip]),
                "active_peak_hz": fpeak,
                "passive_peak_hz": float(freqs[ip]),
                "nonzero_active_peak": False,
            })
            continue

        tone_f = np.asarray([0.35 * fpeak, 0.60 * fpeak, fpeak, 2.5 * fpeak, 6.0 * fpeak], dtype=float)
        ga_soma, gp_soma, ga_prox, gp_prox = [], [], [], []
        for f in tone_f:
            Xa = model.solve_inputs(float(f), [tip], active=True)[:, 0]
            Xp = model.solve_inputs(float(f), [tip], active=False)[:, 0]
            ga_soma.append(abs(Xa[0]))
            gp_soma.append(abs(Xp[0]))
            ga_prox.append(abs(Xa[prox]))
            gp_prox.append(abs(Xp[prox]))

        ga_soma = np.asarray(ga_soma)
        gp_soma = np.asarray(gp_soma)
        ga_prox = np.asarray(ga_prox)
        gp_prox = np.asarray(gp_prox)
        pa_soma = target_fraction(ga_soma)
        pp_soma = target_fraction(gp_soma)
        pa_prox = target_fraction(ga_prox)
        pp_prox = target_fraction(gp_prox)

        rows.append({
            "tip": int(tip),
            "proximal_branch_node": int(prox),
            "section_type": int(tree.section_types[tip]),
            "path_distance_um": float(model.path_distance_um[tip]),
            "active_peak_hz": fpeak,
            "active_peak_over_dc_soma": float(active_soma[ia, j] / max(active_soma[0, j], 1e-300)),
            "passive_peak_hz": float(freqs[ip]),
            "passive_peak_over_dc_soma": float(passive_soma[ip, j] / max(passive_soma[0, j], 1e-300)),
            "nonzero_active_peak": True,
            "tones_hz": [float(x) for x in tone_f],
            "soma": {
                "active_target_fraction": pa_soma,
                "passive_target_fraction": pp_soma,
                "active_over_passive_purity": float(pa_soma / max(pp_soma, 1e-15)),
                "active_total_gain_power": float(np.sum(ga_soma**2)),
                "passive_total_gain_power": float(np.sum(gp_soma**2)),
            },
            "proximal_branch": {
                "active_target_fraction": pa_prox,
                "passive_target_fraction": pp_prox,
                "active_over_passive_purity": float(pa_prox / max(pp_prox, 1e-15)),
                "active_total_gain_power": float(np.sum(ga_prox**2)),
                "passive_total_gain_power": float(np.sum(gp_prox**2)),
            },
        })

    valid = [r for r in rows if r.get("nonzero_active_peak")]
    summary = {
        "nodes": int(model.n),
        "selected_tips": len(tips),
        "nonzero_active_peaks": len(valid),
        "passive_dc_peaks": int(sum(r.get("passive_peak_hz") == 0.0 for r in rows)),
        "active_peak_hz_min": float(min((r["active_peak_hz"] for r in valid), default=0.0)),
        "active_peak_hz_max": float(max((r["active_peak_hz"] for r in valid), default=0.0)),
        "median_active_peak_over_dc": float(np.median([r["active_peak_over_dc_soma"] for r in valid])) if valid else 0.0,
        "median_soma_active_target_fraction": float(np.median([r["soma"]["active_target_fraction"] for r in valid])) if valid else 0.0,
        "median_soma_passive_target_fraction": float(np.median([r["soma"]["passive_target_fraction"] for r in valid])) if valid else 0.0,
        "median_proximal_active_target_fraction": float(np.median([r["proximal_branch"]["active_target_fraction"] for r in valid])) if valid else 0.0,
        "median_proximal_passive_target_fraction": float(np.median([r["proximal_branch"]["passive_target_fraction"] for r in valid])) if valid else 0.0,
    }

    receipt = {
        "gate": 10,
        "question": "Does a fixed quasi-active purifier remain nontrivial on several paths of the real Operaattori human L2/3 morphology?",
        "source": {
            "repository": "ido4848/FCI",
            "commit": SOURCE_COMMIT,
            "path": SOURCE_REL,
            "same_pinned_morphology_as": "anttiluode/Operaattori Gate 10",
        },
        "membrane": {
            "Ra_ohm_cm": model.ra,
            "Cm_F_per_cm2": model.cm,
            "g_leak_S_per_cm2": model.gm,
            "quasi_active_g_over_leak": model.q_over_leak,
            "tau_q_s": model.tau,
            "axon_quasi_active": False,
            "soma_quasi_active": False,
            "warning": "phenomenological restorative density; not fitted HCN expression",
        },
        "summary": summary,
        "paths": rows,
        "claim_boundary": (
            "This gate tests whether real morphology destroys or preserves the operator-level purifier mechanism under one frozen quasi-active parameterization. "
            "It does not show that the real cell has these channel densities or that biological dendrites achieve the simulated purity."
        ),
    }
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")

    print(json.dumps(receipt, indent=2))
    assert model.n > 1000
    assert len(tips) >= 4
    assert np.all(np.isfinite(active_soma))
    assert np.all(np.isfinite(passive_soma))
    return receipt


if __name__ == "__main__":
    run()
