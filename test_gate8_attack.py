from gate8_attack import run


def test_gate7_is_not_one_threshold_or_one_scale_point():
    r = run()
    g8 = r["G8_threshold_attack"]
    scale = r["G8_scale_attack"]

    assert g8["robust_threshold_count"] >= 3
    assert g8["frozen_threshold_is_in_robust_window"]
    assert g8["robust_threshold_min"] <= 0.0115 + 1e-12
    assert g8["robust_threshold_max"] >= 0.012
    assert scale["minimum_route_selectivity_ratio"] > 8.0
