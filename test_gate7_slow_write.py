from gate7_slow_write import run


def test_persistent_local_write_survives_fast_state_erasure():
    r = run()
    g7 = r["G7_persistent_local_write"]
    recall = r["G7_after_fast_state_erasure"]

    assert g7["individual_A_write"] == 0.0
    assert g7["individual_B_write"] == 0.0
    assert g7["matched_over_separated_write"] > 10.0
    assert recall["relative_B_change_matched"] > 0.05
    assert recall["relative_B_change_separated"] < 0.02
    assert recall["branch_vs_soma_relative_visibility"] > 5.0
