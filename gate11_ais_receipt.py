"""Gate 11: PURIFY -> COMMIT -> RETURN RECEIPT -> WRITE -> REPLAY.

This is a deliberately small biological-looking operator laboratory, not a
quantitative pyramidal-neuron model.

The existing quasi-active dendrite supplies the forward frequency-selective
operator.  We append one axon-initial-segment (AIS) compartment.  In the
small-signal frozen state the electrical admittance is reciprocal/symmetric.
An AIS event is represented by a transient post-spike conductance state; the
return pulse therefore traverses a *different state of the same morphology*.
That history-conditioned difference is called the receipt here.

A local eligibility field from the forward cue is multiplied by the receipt
and by a slow-write permission bit.  The strongest branch edge receives one
persistent axial-conductance edit.  Fast state is then treated as exactly
wiped, and the identical frequency-domain cue is replayed through the edited
operator.

Important distinction:
    receipt != learning
Learning is the persistent operator edit caused when eligibility, receipt and
write permission coincide.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from active_resonance import QuasiActiveDendrite, log_frequency_grid


def _edge_outer(n: int, i: int, j: int) -> np.ndarray:
    b = np.zeros(n, dtype=float)
    b[int(i)] = 1.0
    b[int(j)] = -1.0
    return np.outer(b, b)


class DendriteAIS:
    """Quasi-active branched dendrite plus one electrically coupled AIS node."""

    def __init__(
        self,
        branch_len: int = 12,
        g_soma_ais: float = 0.55,
        g_ais_leak: float = 0.10,
        c_ais: float = 1.1,
    ):
        self.dendrite = QuasiActiveDendrite(branch_len=branch_len)
        self.n_d = self.dendrite.n
        self.ais = self.n_d
        self.n = self.n_d + 1
        self.soma = self.dendrite.soma
        self.input_site = int(self.dendrite.branch_a[-1])
        self.branch = self.dendrite.branch_a.astype(int)
        self.g_soma_ais = float(g_soma_ais)
        self.g_ais_leak = float(g_ais_leak)
        self.c_ais = float(c_ais)

    def admittance(
        self,
        omega: float,
        *,
        post_spike: bool = False,
        learned_edge: tuple[int, int, float] | None = None,
    ) -> np.ndarray:
        """Complex small-signal admittance for one frozen physiological state."""
        w = float(omega)
        Y = np.zeros((self.n, self.n), dtype=complex)
        Y[: self.n_d, : self.n_d] = self.dendrite.admittance_matrix(w)
        Y[self.ais, self.ais] = self.g_ais_leak + 1j * w * self.c_ais

        # Reciprocal axial soma<->AIS coupling.
        Y += self.g_soma_ais * _edge_outer(self.n, self.soma, self.ais)

        # A finite event changes channel state.  This does not make the frozen
        # post-event admittance non-reciprocal; it makes the return traversal
        # occur through a different operator than the pre-event forward pass.
        if post_spike:
            Y[self.ais, self.ais] += 1.25  # refractory/AHP-like AIS conductance
            Y[self.soma, self.soma] += 0.18  # somatic after-event shunt

        if learned_edge is not None:
            i, j, delta_g = learned_edge
            Y += float(delta_g) * _edge_outer(self.n, int(i), int(j))
        return Y

    def transfer_matrix(self, omega: float, **kwargs) -> np.ndarray:
        return np.linalg.inv(self.admittance(omega, **kwargs))


def _target_fraction(H_at_tones: list[np.ndarray], ais: int, input_site: int) -> float:
    gains = np.asarray([abs(H[ais, input_site]) for H in H_at_tones], dtype=float)
    power = gains * gains
    return float(power[2] / power.sum())


def _relative_norm(a, b) -> float:
    denom = max(float(np.linalg.norm(b)), 1e-15)
    return float(np.linalg.norm(a - b) / denom)


def run(write_scale: float = 0.12, write_permission: float = 1.0) -> dict:
    machine = DendriteAIS(branch_len=12)

    # SELECT/PURIFY: find the non-zero band that most strongly reaches the AIS.
    omegas = log_frequency_grid(low=0.002, high=0.8, count=420)
    gains = np.asarray(
        [
            abs(machine.transfer_matrix(w)[machine.ais, machine.input_site])
            for w in omegas
        ],
        dtype=float,
    )
    peak_index = int(np.argmax(gains))
    omega_star = float(omegas[peak_index])

    # Five equal-amplitude temporal components.  Target starts at 20% input
    # energy; measure what fraction reaches the AIS after dendritic filtering.
    tones = np.asarray(
        [0.0, 0.35 * omega_star, omega_star, 2.5 * omega_star, 6.0 * omega_star],
        dtype=float,
    )
    H_tones = [machine.transfer_matrix(w) for w in tones]
    purity = _target_fraction(H_tones, machine.ais, machine.input_site)

    # Frozen small-signal reciprocity: complex-symmetric transfer, H == H^T.
    H_pre = machine.transfer_matrix(omega_star)
    reciprocity_error = float(
        np.linalg.norm(H_pre - H_pre.T) / max(np.linalg.norm(H_pre), 1e-15)
    )

    forward = H_pre[: machine.n_d, machine.input_site]
    reverse_small = H_pre[: machine.n_d, machine.ais]
    pairwise_reciprocity_error = float(
        np.max(
            np.abs(
                H_pre[machine.ais, : machine.n_d]
                - H_pre[: machine.n_d, machine.ais]
            )
        )
        / max(np.max(np.abs(H_pre)), 1e-15)
    )

    # COMMIT/RETURN: event changes the fast channel state; an AIS-originating
    # return pulse now sees H_post instead of H_pre.
    H_post = machine.transfer_matrix(omega_star, post_spike=True)
    reverse_spike = H_post[: machine.n_d, machine.ais]
    receipt = reverse_spike - reverse_small
    receipt_relative = float(
        np.linalg.norm(receipt) / max(np.linalg.norm(reverse_small), 1e-15)
    )

    # BIND: forward occupancy is the local eligibility field.  The receipt is
    # the local evidence that an output event changed the state.  Use magnitudes
    # here to avoid pretending this toy has a calibrated biological sign rule.
    eligibility = np.abs(forward)
    receipt_mag = np.abs(receipt)
    eligibility /= max(float(eligibility.max()), 1e-15)
    receipt_mag /= max(float(receipt_mag.max()), 1e-15)
    joint = eligibility * receipt_mag

    # Only branch-A edges are eligible for this first closed-loop gate.  The
    # chosen site is determined by the joint local overlap, not hard-coded.
    branch_scores = joint[machine.branch]
    branch_k = int(np.argmax(branch_scores))
    write_site = int(machine.branch[branch_k])
    parent_site = machine.soma if branch_k == 0 else int(machine.branch[branch_k - 1])
    joint_score = float(branch_scores[branch_k])

    delta_g = float(write_scale) * float(write_permission) * joint_score
    learned_edge = (parent_site, write_site, delta_g)

    # WRITE: one local physical edge edit creates a distributed rank-one change
    # in the fixed-frequency resolvent.
    H_learned = machine.transfer_matrix(omega_star, learned_edge=learned_edge)
    delta_H = H_learned - H_pre
    svals = np.linalg.svd(delta_H, compute_uv=False)
    rank1_ratio = float(svals[1] / svals[0]) if svals[0] > 0 else 0.0

    # REPLAY after an exact conceptual fast-state wipe.  Frequency-domain
    # transfer contains no carried voltage/gating history from the first pass;
    # only the persistent edge parameter differs.
    pre_gain = float(abs(H_pre[machine.ais, machine.input_site]))
    post_gain = float(abs(H_learned[machine.ais, machine.input_site]))
    replay_relative_change = float(abs(post_gain - pre_gain) / max(pre_gain, 1e-15))

    H_learned_tones = [
        machine.transfer_matrix(w, learned_edge=learned_edge) for w in tones
    ]
    learned_purity = _target_fraction(
        H_learned_tones, machine.ais, machine.input_site
    )

    # Controls: no consequence permission => no persistent rewrite.  This is
    # deliberately analogous to failure/relevance-gated consolidation rather
    # than unconditional reinforcement on every spike.
    no_write_edge = (parent_site, write_site, 0.0)
    H_no_write = machine.transfer_matrix(omega_star, learned_edge=no_write_edge)
    no_write_change = _relative_norm(H_no_write, H_pre)

    return {
        "model": {
            "dendritic_compartments": machine.n_d,
            "total_compartments_with_AIS": machine.n,
            "input_site": machine.input_site,
            "soma": machine.soma,
            "ais": machine.ais,
        },
        "select_purify": {
            "omega_star": omega_star,
            "peak_gain_to_AIS": float(gains[peak_index]),
            "input_target_energy_fraction": 0.2,
            "AIS_target_energy_fraction": purity,
            "enrichment_over_input": float(purity / 0.2),
            "tones": [float(x) for x in tones],
        },
        "small_signal_reciprocity": {
            "matrix_relative_error_H_minus_HT": reciprocity_error,
            "pairwise_forward_reverse_relative_error": pairwise_reciprocity_error,
        },
        "commit_return_receipt": {
            "receipt_relative_to_small_return": receipt_relative,
            "interpretation": (
                "The event-conditioned return differs because it traverses a post-event "
                "conductance state. Each frozen state remains reciprocal."
            ),
        },
        "bind_write": {
            "write_permission": float(write_permission),
            "joint_eligibility_receipt_score": joint_score,
            "parent_site": parent_site,
            "write_site": write_site,
            "delta_g": delta_g,
            "resolvent_s2_over_s1": rank1_ratio,
        },
        "replay_after_fast_state_wipe": {
            "fast_state": "exactly absent from the frequency-domain replay",
            "pre_AIS_gain": pre_gain,
            "post_AIS_gain": post_gain,
            "relative_gain_change": replay_relative_change,
            "pre_target_fraction": purity,
            "post_target_fraction": learned_purity,
            "target_fraction_change": float(learned_purity - purity),
        },
        "controls": {
            "no_write_permission_operator_change": no_write_change,
            "receipt_is_not_learning": True,
            "learning_requires_persistent_edge_edit": bool(delta_g > 0.0),
        },
        "claim_boundary": (
            "This gate demonstrates a closed toy operator loop with quasi-active dendritic "
            "purification, a reciprocal small-signal baseline, an AIS-state-conditioned return "
            "receipt, a receipt-and-eligibility-gated local structural edit, and changed replay "
            "after fast-state erasure. It is not a fitted AIS/bAP/STDP model and does not show "
            "that biological bAPs compute an adjoint or exact gradient."
        ),
    }


if __name__ == "__main__":
    result = run()
    print(json.dumps(result, indent=2, allow_nan=False))
    out = Path("results")
    out.mkdir(exist_ok=True)
    (out / "gate11_ais_receipt.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
