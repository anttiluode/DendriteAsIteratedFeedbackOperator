from nonlinear_experiment import run


def test_gate6_relationship_verdict_passes():
    r = run()
    assert r["verdict"]["passes"]
    assert r["summary"]["matched_over_largest_control"] > 1.0


def test_gate6_is_not_single_input_thresholding():
    r = run()
    knee = r["model"]["regenerative_knee"]
    assert r["cases"]["matched"]["linear_hotspot_peak"] > knee
    assert r["cases"]["single_A"]["linear_hotspot_peak"] < knee
    assert r["cases"]["single_B"]["linear_hotspot_peak"] < knee
