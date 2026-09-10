import numpy as np

from active_resonance import log_frequency_grid
from gate11_ais_receipt import DendriteAIS
from gate12_natural_length import LengthOperator, tone_metrics, grounded_score


def _fixture():
    machine = DendriteAIS(branch_len=12)
    parent = int(machine.branch[-2])
    site = int(machine.branch[-1])
    op = LengthOperator(machine, parent, site)
    omegas = log_frequency_grid(low=0.002, high=0.8, count=420)
    gains = np.asarray(
        [abs(machine.transfer_matrix(w)[machine.ais, machine.input_site]) for w in omegas]
    )
    omega_star = float(omegas[int(np.argmax(gains))])
    tones = np.asarray([0.0, 0.35 * omega_star, omega_star, 2.5 * omega_star, 6.0 * omega_star])
    return machine, op, tones


def test_unit_length_is_exact_baseline_operator():
    machine, op, tones = _fixture()
    for w in tones:
        a = op.admittance(float(w), 1.0)
        b = machine.admittance(float(w))
        assert np.linalg.norm(a - b) < 1e-12


def test_purity_can_improve_while_signal_collapses():
    _, op, tones = _fixture()
    base = tone_metrics(op, tones, 1.0)
    isolated = tone_metrics(op, tones, 50.0)
    assert isolated["target_fraction"] > base["target_fraction"]
    assert isolated["target_gain"] < 0.01 * base["target_gain"]


def test_grounded_observer_rejects_near_isolation():
    _, op, tones = _fixture()
    base = tone_metrics(op, tones, 1.0)
    noise = 0.01 * base["target_power"]
    score_base = grounded_score(base, noise)
    score_isolated = grounded_score(tone_metrics(op, tones, 50.0), noise)
    assert score_base > score_isolated
