"""Gate 6 receipt: does branch-selected resonance control a local nonlinear event?"""

from __future__ import annotations

import json
from pathlib import Path

from nonlinear_sieve import ResonantNonlinearSieve


def run():
    sieve = ResonantNonlinearSieve()
    cases = sieve.run_all()

    r = {k: v["nonlinear_residual_soma_rms"] for k, v in cases.items()}
    matched = r["matched"]
    offband_max = max(r["wrong_low_carrier"], r["wrong_high_carrier"])
    single_max = max(r["single_A"], r["single_B"])
    control_max = max(
        r["quadrature"],
        r["antiphase"],
        offband_max,
        r["spatially_separated"],
        single_max,
    )

    verdict = {
        "matched_linear_hotspot_crosses_knee": cases["matched"]["linear_hotspot_peak"] > sieve.threshold,
        "single_A_stays_below_knee": cases["single_A"]["linear_hotspot_peak"] < sieve.threshold,
        "single_B_stays_below_knee": cases["single_B"]["linear_hotspot_peak"] < sieve.threshold,
        "matched_beats_quadrature": matched > r["quadrature"],
        "matched_beats_antiphase": matched > r["antiphase"],
        "matched_beats_offband": matched > offband_max,
        "matched_beats_spatial_separation": matched > r["spatially_separated"],
        "matched_beats_singles": matched > single_max,
    }
    verdict["passes"] = bool(all(verdict.values()))

    return {
        "model": {
            "target_branch_peak_omega": sieve.omega,
            "hotspot": sieve.hotspot,
            "same_branch_sites": [sieve.site_a, sieve.site_b],
            "separated_control_site": sieve.separated_site_b,
            "relative_input_amplitudes": [sieve.rel_amp_a, sieve.rel_amp_b],
            "drive_scale": sieve.drive_scale,
            "regenerative_knee": sieve.threshold,
            "regenerative_slope": sieve.slope,
            "regenerative_Imax": sieve.i_max,
            "note": "bounded generic regenerative hotspot; not an NMDA kinetic fit",
        },
        "cases": cases,
        "summary": {
            "matched_residual_soma_rms": matched,
            "largest_control_residual_soma_rms": control_max,
            "matched_over_largest_control": matched / max(control_max, 1e-15),
            "matched_over_quadrature": matched / max(r["quadrature"], 1e-15),
            "matched_over_antiphase": matched / max(r["antiphase"], 1e-15),
            "matched_over_offband": matched / max(offband_max, 1e-15),
            "matched_over_spatial_separation": matched / max(r["spatially_separated"], 1e-15),
            "matched_over_largest_single": matched / max(single_max, 1e-15),
        },
        "verdict": verdict,
    }


if __name__ == "__main__":
    receipt = run()
    print(json.dumps(receipt, indent=2))
    out = Path("results")
    out.mkdir(exist_ok=True)
    (out / "nonlinear_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if not receipt["verdict"]["passes"]:
        raise SystemExit("Gate 6 nonlinear resonant-sieve verdict failed")
