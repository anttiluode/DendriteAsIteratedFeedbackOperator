from gate11_ais_receipt import run


def test_gate11_closes_operator_loop():
    r = run()

    assert r["select_purify"]["omega_star"] > 0.0
    assert r["select_purify"]["AIS_target_energy_fraction"] > 0.20

    assert (
        r["small_signal_reciprocity"]["matrix_relative_error_H_minus_HT"]
        < 1e-10
    )
    assert (
        r["small_signal_reciprocity"]["pairwise_forward_reverse_relative_error"]
        < 1e-10
    )

    assert r["commit_return_receipt"]["receipt_relative_to_small_return"] > 1e-3

    assert r["bind_write"]["joint_eligibility_receipt_score"] > 0.0
    assert r["bind_write"]["delta_g"] > 0.0
    assert r["bind_write"]["resolvent_s2_over_s1"] < 1e-10

    assert r["replay_after_fast_state_wipe"]["relative_gain_change"] > 1e-6
    assert r["controls"]["no_write_permission_operator_change"] < 1e-12
    assert r["controls"]["receipt_is_not_learning"] is True


def test_gate11_permission_can_block_persistent_write():
    r = run(write_permission=0.0)
    assert r["bind_write"]["delta_g"] == 0.0
    assert r["replay_after_fast_state_wipe"]["relative_gain_change"] < 1e-12
    assert r["controls"]["learning_requires_persistent_edge_edit"] is False
