from gate9_causal_write import run


def test_threshold_free_causal_component_survives_fast_state_erasure():
    r = run()
    write = r["G9_threshold_free_write"]
    recall = r["G9_after_fast_state_erasure"]
    sweep = r["G9_eta_attack"]

    assert write["individual_A_write"] > 0.0
    assert write["individual_B_write"] > 0.0
    assert write["matched_over_separated_collision_write"] > 1.5
    assert recall["matched_over_separated_B_chi_norm"] > 1.5
    assert recall["matched_over_separated_B_chi_peak"] > 1.5
    assert sweep["minimum_B_chi_selectivity"] > 1.5
