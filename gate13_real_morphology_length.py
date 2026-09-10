"""Gate 13: does the grounded natural-length idea survive a real dendritic arbor?

This gate transfers the Gate-12 structural question onto a pinned human L2/3
pyramidal morphology used by Operaattori.  The morphology is real; the membrane
model remains a deliberately simple, uniform quasi-active cable model.

We ask a narrow question:

    If one real dendritic segment is allowed to change length while its radius
    and membrane densities stay fixed, does a bounded-observer objective have a
    finite local fixed point?

The segment is not chosen by looking for the best length derivative.  It is
chosen by a forward/reverse overlap at the baseline resonant frequency:

    |forward occupancy| x |reciprocal return|

along the actual soma-to-input path.  This keeps the Gate-11 intuition that the
same geometry which carries influence toward the output also identifies where
an output-originating receipt can return strongly.

Important claim boundary:
- real morphology: yes (pinned ASC, cell 1125)
- quantitative ion-channel fit: no
- biological growth rule: no
- literal acoustic wavelength matching: no
- exact bAP / adjoint: no

The gate is about geometry compiling a frequency-dependent operator and about
whether a local length update can stop at a useful, observer-grounded geometry.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import urllib.request

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


SOURCE_COMMIT = "75ad8b4d81a7f51bf888b30650c543592340db06"
SOURCE_NAME = "2013_03_06_cell11_1125_H41_06.asc"
SOURCE_REL = (
    "simulating_neurons/neuron_models/human/eyal/"
    "Human_L23_PC_0603_11_937_Eyal_passive_dends_simple_soma/"
    "morphologies/" + SOURCE_NAME
)
SOURCE_URL = "https://raw.githubusercontent.com/ido4848/FCI/" + SOURCE_COMMIT + "/" + SOURCE_REL

# Simple, explicit membrane/cable constants.  Geometry is in micrometres,
# converted below to centimetres for these conventional cable quantities.
RA_OHM_CM = 150.0
CM_F_PER_CM2 = 1.0e-6
GLEAK_S_PER_CM2 = 1.0e-4
GQ_S_PER_CM2 = 1.5e-3
TAU_Q_S = 0.050
G_SOMA_AIS_S = 2.0e-7
G_AIS_LEAK_S = 2.0e-9
C_AIS_F = 2.0e-11


@dataclass(frozen=True)
class PointTree:
    positions_um: np.ndarray
    parents: np.ndarray
    radii_um: np.ndarray
    section_types: np.ndarray


def _enum_int(value: object) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        raw = getattr(value, "value", None)
        return -1 if raw is None else int(raw)


def download_source(dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 100_000:
        return dest
    req = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "DendriteGate13/1.0"})
    with urllib.request.urlopen(req, timeout=60) as response, dest.open("wb") as f:
        f.write(response.read())
    return dest


def load_dendritic_tree(path: Path, duplicate_tol_um: float = 1e-5) -> PointTree:
    """Load the pinned morphology and retain soma + dendritic section types 3/4."""
    from morphio import Morphology

    morph = Morphology(str(path))
    soma_points = np.asarray(morph.soma.points, dtype=float)
    soma_diameters = np.asarray(morph.soma.diameters, dtype=float)
    if len(soma_points):
        root = np.mean(soma_points, axis=0)
        root_radius = float(np.mean(soma_diameters) / 2.0) if len(soma_diameters) else 5.0
    else:
        root = np.zeros(3, dtype=float)
        root_radius = 5.0

    positions = [np.asarray(root, dtype=float)]
    parents = [-1]
    radii = [max(root_radius, 0.05)]
    section_types = [1]
    section_end: dict[int, int] = {}
    skipped_sections: set[int] = set()

    for section in morph.iter():
        sid = int(section.id)
        stype = _enum_int(section.type)
        if stype not in (3, 4):
            skipped_sections.add(sid)
            continue

        if section.is_root:
            current_parent = 0
        else:
            pid = int(section.parent.id)
            if pid in skipped_sections:
                # Dendrites in this morphology should connect to soma or another
                # dendritic section.  Fail loudly rather than silently inventing
                # a topological shortcut.
                raise ValueError(f"dendritic section {sid} has skipped parent section {pid}")
            if pid not in section_end:
                raise ValueError(f"parent section {pid} not available before child {sid}")
            current_parent = section_end[pid]

        pts = np.asarray(section.points, dtype=float)
        diams = np.asarray(section.diameters, dtype=float)
        if len(pts) != len(diams):
            raise ValueError("point/diameter mismatch")
        start = 0
        if len(pts) and np.linalg.norm(pts[0] - positions[current_parent]) <= duplicate_tol_um:
            start = 1
        for j in range(start, len(pts)):
            p = np.asarray(pts[j], dtype=float)
            if np.linalg.norm(p - positions[current_parent]) <= duplicate_tol_um:
                continue
            idx = len(positions)
            positions.append(p)
            parents.append(current_parent)
            radii.append(max(float(diams[j]) / 2.0, 0.05))
            section_types.append(stype)
            current_parent = idx
        section_end[sid] = current_parent

    tree = PointTree(
        positions_um=np.asarray(positions, dtype=float),
        parents=np.asarray(parents, dtype=np.int64),
        radii_um=np.asarray(radii, dtype=float),
        section_types=np.asarray(section_types, dtype=np.int64),
    )
    if len(tree.positions_um) < 1000:
        raise ValueError("unexpectedly small dendritic tree")
    return tree


def edge_lengths_um(tree: PointTree) -> np.ndarray:
    out = np.zeros(len(tree.parents), dtype=float)
    for i in range(1, len(out)):
        out[i] = np.linalg.norm(tree.positions_um[i] - tree.positions_um[int(tree.parents[i])])
    return out


def child_counts(parents: np.ndarray) -> np.ndarray:
    out = np.zeros(len(parents), dtype=np.int64)
    for i in range(1, len(parents)):
        out[int(parents[i])] += 1
    return out


def choose_distal_tip(tree: PointTree, lengths_um: np.ndarray) -> tuple[int, np.ndarray, float]:
    """Choose the dendritic tip with the greatest cable distance from soma."""
    dist = np.zeros(len(tree.parents), dtype=float)
    for i in range(1, len(dist)):
        dist[i] = dist[int(tree.parents[i])] + lengths_um[i]
    tips = np.flatnonzero(child_counts(tree.parents) == 0)
    tips = tips[tips != 0]
    tip = int(tips[int(np.argmax(dist[tips]))])
    path = []
    j = tip
    while j != 0:
        path.append(j)
        j = int(tree.parents[j])
    path.append(0)
    path = np.asarray(path[::-1], dtype=np.int64)
    return tip, path, float(dist[tip])


def _membrane_area_cm2(radius_um: float, length_um: float) -> float:
    # Cylinder lateral area in um^2, then um^2 -> cm^2.
    return float(2.0 * np.pi * radius_um * length_um * 1.0e-8)


def build_operator(
    tree: PointTree,
    lengths_um: np.ndarray,
    *,
    edit_node: int | None = None,
    length_scale: float = 1.0,
) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray, int]:
    """Compile real geometry into a simple quasi-active cable + AIS operator."""
    if length_scale <= 0:
        raise ValueError("length_scale must be positive")
    n_d = len(tree.parents)
    ais = n_d
    n = n_d + 1

    lengths = np.asarray(lengths_um, dtype=float).copy()
    if edit_node is not None:
        if int(edit_node) <= 0 or int(edit_node) >= n_d:
            raise ValueError("edit_node must be a non-root dendritic node")
        lengths[int(edit_node)] *= float(length_scale)

    leak = np.zeros(n, dtype=float)
    cap = np.zeros(n, dtype=float)
    gq = np.zeros(n, dtype=float)

    # Soma root: spherical approximation, only to close the cable electrically.
    soma_area_cm2 = 4.0 * np.pi * tree.radii_um[0] ** 2 * 1.0e-8
    leak[0] = GLEAK_S_PER_CM2 * soma_area_cm2
    cap[0] = CM_F_PER_CM2 * soma_area_cm2
    gq[0] = GQ_S_PER_CM2 * soma_area_cm2

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []

    for i in range(1, n_d):
        p = int(tree.parents[i])
        L_um = max(float(lengths[i]), 1e-6)
        r_um = max(float(tree.radii_um[i]), 0.05)
        area = _membrane_area_cm2(r_um, L_um)
        leak[i] = GLEAK_S_PER_CM2 * area
        cap[i] = CM_F_PER_CM2 * area
        gq[i] = GQ_S_PER_CM2 * area

        # Axial conductance through a cylindrical segment.  Use mean endpoint
        # radius to avoid treating every taper as a child-radius discontinuity.
        r_edge_um = max(0.5 * (tree.radii_um[p] + tree.radii_um[i]), 0.05)
        r_cm = r_edge_um * 1.0e-4
        L_cm = L_um * 1.0e-4
        g_ax = float(np.pi * r_cm * r_cm / (RA_OHM_CM * L_cm))
        rows.extend([i, p, i, p])
        cols.extend([i, p, p, i])
        vals.extend([g_ax, g_ax, -g_ax, -g_ax])

    # AIS-like bounded output node.
    leak[ais] = G_AIS_LEAK_S
    cap[ais] = C_AIS_F
    rows.extend([0, ais, 0, ais])
    cols.extend([0, ais, ais, 0])
    vals.extend([G_SOMA_AIS_S, G_SOMA_AIS_S, -G_SOMA_AIS_S, -G_SOMA_AIS_S])

    # Add membrane leak on the diagonal.
    for i in range(n):
        rows.append(i)
        cols.append(i)
        vals.append(float(leak[i]))

    G = sparse.coo_matrix((vals, (rows, cols)), shape=(n, n), dtype=float).tocsr()
    G.sum_duplicates()
    return G, cap, gq, ais


def transfer_vector(
    tree: PointTree,
    lengths_um: np.ndarray,
    omega_rad_s: float,
    source_node: int,
    *,
    edit_node: int | None = None,
    length_scale: float = 1.0,
) -> tuple[np.ndarray, int]:
    G, cap, gq, ais = build_operator(
        tree, lengths_um, edit_node=edit_node, length_scale=length_scale
    )
    w = float(omega_rad_s)
    yq = gq / (1.0 + 1j * w * TAU_Q_S)
    A = G.astype(complex) + sparse.diags(1j * w * cap + yq, format="csr")
    b = np.zeros(A.shape[0], dtype=complex)
    b[int(source_node)] = 1.0
    return np.asarray(spsolve(A.tocsc(), b), dtype=complex), ais


def tone_transfer(
    tree: PointTree,
    lengths_um: np.ndarray,
    tones: np.ndarray,
    input_node: int,
    *,
    edit_node: int | None = None,
    length_scale: float = 1.0,
) -> np.ndarray:
    out = []
    for w in tones:
        v, ais = transfer_vector(
            tree,
            lengths_um,
            float(w),
            input_node,
            edit_node=edit_node,
            length_scale=length_scale,
        )
        out.append(v[ais])
    return np.asarray(out, dtype=complex)


def purity(z: np.ndarray) -> float:
    p = np.abs(np.asarray(z, dtype=complex)) ** 2
    return float(p[2] / max(float(np.sum(p)), 1e-300))


def grounded_score(z: np.ndarray, noise_power: float) -> float:
    p = np.abs(np.asarray(z, dtype=complex)) ** 2
    target = float(p[2])
    off = float(np.sum(p) - target)
    return target / max(off + float(noise_power), 1e-300)


def discrete_local_fixed_point(values: np.ndarray, start_index: int) -> tuple[int, list[int]]:
    """Move to the better immediate grid neighbour until neither neighbour wins."""
    values = np.asarray(values, dtype=float)
    i = int(start_index)
    history = [i]
    for _ in range(len(values) + 2):
        candidates = [i]
        if i > 0:
            candidates.append(i - 1)
        if i + 1 < len(values):
            candidates.append(i + 1)
        j = max(candidates, key=lambda k: values[k])
        if j == i:
            return i, history
        i = int(j)
        history.append(i)
    raise RuntimeError("local fixed-point walk failed to terminate")


def run() -> dict:
    out_dir = Path("results")
    source = download_source(out_dir / "_source" / SOURCE_NAME)
    tree = load_dendritic_tree(source)
    lengths = edge_lengths_um(tree)
    tip, path, path_length = choose_distal_tip(tree, lengths)

    # Baseline frequency sweep.  The target carrier is fixed once from the
    # unedited real morphology and is not redefined after growth.
    omegas = np.geomspace(0.5, 500.0, 120)
    gains = []
    for w in omegas:
        v, ais = transfer_vector(tree, lengths, float(w), tip)
        gains.append(abs(v[ais]))
    gains = np.asarray(gains, dtype=float)
    peak_i = int(np.argmax(gains))
    omega_star = float(omegas[peak_i])
    tones = np.asarray([0.0, 0.35 * omega_star, omega_star, 2.5 * omega_star, 6.0 * omega_star])
    z0 = tone_transfer(tree, lengths, tones, tip)
    p0 = np.abs(z0) ** 2

    # Gate-11-inspired WHERE signal: baseline forward occupancy from the distal
    # input multiplied by reciprocal return from the AIS, restricted to the
    # actual soma-to-input path.  Exclude root and terminal tip so the edit is
    # an internal path segment rather than a trivial endpoint isolation knob.
    forward, ais = transfer_vector(tree, lengths, omega_star, tip)
    reverse, _ = transfer_vector(tree, lengths, omega_star, ais)
    joint = np.abs(forward[: len(tree.parents)]) * np.abs(reverse[: len(tree.parents)])
    joint /= max(float(np.max(joint[path])), 1e-300)
    eligible_path = path[1:-1]
    edit_node = int(eligible_path[int(np.argmax(joint[eligible_path]))])
    edit_parent = int(tree.parents[edit_node])

    # One common geometry sweep is enough to evaluate many observer noise floors.
    scales = np.geomspace(0.2, 5.0, 101)
    transfer_grid = np.asarray(
        [tone_transfer(tree, lengths, tones, tip, edit_node=edit_node, length_scale=float(s)) for s in scales],
        dtype=complex,
    )
    powers = np.abs(transfer_grid) ** 2
    purities = powers[:, 2] / np.maximum(np.sum(powers, axis=1), 1e-300)
    start_i = int(np.argmin(np.abs(np.log(scales))))

    purity_fp_i, purity_history = discrete_local_fixed_point(purities, start_i)

    noise_fractions = [0.001, 0.003, 0.01, 0.03, 0.1]
    rows = []
    reference_i = None
    for frac in noise_fractions:
        noise = float(frac * p0[2])
        scores = powers[:, 2] / np.maximum(np.sum(powers, axis=1) - powers[:, 2] + noise, 1e-300)
        fp_i, history = discrete_local_fixed_point(scores, start_i)
        broad_i = int(np.argmax(scores))
        row = {
            "observer_noise_fraction_of_baseline_target_power": frac,
            "noise_floor_power": noise,
            "fixed_point_scale": float(scales[fp_i]),
            "fixed_point_length_um": float(lengths[edit_node] * scales[fp_i]),
            "fixed_point_score": float(scores[fp_i]),
            "fixed_point_is_interior": bool(0 < fp_i < len(scales) - 1),
            "broad_best_scale": float(scales[broad_i]),
            "broad_best_is_interior": bool(0 < broad_i < len(scales) - 1),
            "target_fraction": float(purities[fp_i]),
            "target_gain_vs_baseline": float(abs(transfer_grid[fp_i, 2]) / max(abs(z0[2]), 1e-300)),
            "walk_steps": len(history) - 1,
        }
        rows.append(row)
        if abs(frac - 0.01) < 1e-12:
            reference_i = fp_i

    assert reference_i is not None
    zr = transfer_grid[int(reference_i)]
    scalar = np.vdot(z0, zr) / max(float(np.vdot(z0, z0).real), 1e-300)
    scalar_residual = float(np.linalg.norm(zr - scalar * z0) / max(np.linalg.norm(zr), 1e-300))
    phase_shift = float(np.angle(zr[2] / z0[2]))

    children = child_counts(tree.parents)
    result = {
        "source": {
            "repository": "ido4848/FCI",
            "commit": SOURCE_COMMIT,
            "path": SOURCE_REL,
            "morphology_identifier": "1125",
            "cell": "human L2/3 pyramidal morphology",
        },
        "real_morphology": {
            "dendritic_nodes_including_soma_root": int(len(tree.parents)),
            "dendritic_edges": int(len(tree.parents) - 1),
            "branch_nodes": int(np.sum(children >= 2)),
            "tips": int(np.sum(children == 0) - 1),
            "total_dendritic_cable_um": float(np.sum(lengths)),
            "selected_input_tip": tip,
            "selected_input_path_nodes": int(len(path)),
            "selected_input_path_length_um": path_length,
        },
        "phenomenological_operator": {
            "Ra_ohm_cm": RA_OHM_CM,
            "Cm_F_per_cm2": CM_F_PER_CM2,
            "gleak_S_per_cm2": GLEAK_S_PER_CM2,
            "gq_S_per_cm2": GQ_S_PER_CM2,
            "tau_q_s": TAU_Q_S,
            "omega_star_rad_s": omega_star,
            "peak_is_interior_to_frequency_sweep": bool(0 < peak_i < len(omegas) - 1),
            "baseline_target_fraction": purity(z0),
            "baseline_target_gain": float(abs(z0[2])),
            "tones_rad_s": [float(x) for x in tones],
        },
        "receipt_selected_geometry": {
            "selection_rule": "max |forward_from_tip| * |reciprocal_return_from_AIS| on the real input path",
            "parent_node": edit_parent,
            "edit_node": edit_node,
            "baseline_segment_length_um": float(lengths[edit_node]),
            "segment_radius_um": float(tree.radii_um[edit_node]),
            "normalized_joint_overlap": float(joint[edit_node]),
            "not_selected_by_length_gradient": True,
        },
        "purity_only": {
            "fixed_point_scale": float(scales[purity_fp_i]),
            "fixed_point_is_interior": bool(0 < purity_fp_i < len(scales) - 1),
            "target_fraction": float(purities[purity_fp_i]),
            "target_gain_vs_baseline": float(abs(transfer_grid[purity_fp_i, 2]) / max(abs(z0[2]), 1e-300)),
            "walk_steps": len(purity_history) - 1,
        },
        "grounded_bounded_observer": {
            "score": "target_power / (off_target_power + observer_noise_floor)",
            "rows": rows,
            "reference_noise_fraction": 0.01,
            "reference_fixed_point_scale": float(scales[int(reference_i)]),
            "reference_fixed_point_length_um": float(lengths[edit_node] * scales[int(reference_i)]),
        },
        "operator_change_at_reference_length": {
            "target_transfer_phase_shift_radians": phase_shift,
            "five_tone_best_scalar_residual": scalar_residual,
            "interpretation": "A real segment-length edit changes the complex frequency response and is not generally equivalent to one scalar weight across the five tones.",
        },
        "sigh_bridge": (
            "The fixed point is defined at the structural level: shorter/current/longer is iterated until the bounded observer has no better neighbouring geometry. "
            "This is the geometry analogue of an invariant-state stopping condition, with an absolute observation floor preventing mode purity from being confused with useful transmission."
        ),
        "claim_boundary": (
            "The morphology and segment metrics are real, but the membrane kinetics are a uniform phenomenological quasi-active model. "
            "A finite fixed point is task- and observer-relative. This does not show that biological dendrites execute this search, use this return signal, or grow to acoustic wavelengths."
        ),
    }
    out_dir.mkdir(exist_ok=True)
    (out_dir / "gate13_real_morphology_length.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    run()
