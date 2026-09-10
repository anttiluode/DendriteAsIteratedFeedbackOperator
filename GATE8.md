# Gate 8 — the first parameter attacker

Gate 7 used an engineered plasticity threshold and write scale. A single successful point is too easy to dismiss, so Gate 8 keeps the dynamics fixed and attacks those two knobs.

## Threshold sweep

Eighteen thresholds from `0.007` to `0.0155` are tested at the original write scale. A point counts as robust only if all four conditions hold at once:

```text
individual A+B contamination <= 1%
matched / separated write     >= 3x
matched later B-route change  >= 5%
separated route change        <= 2%
```

Only three sampled thresholds survive:

```text
0.0115
0.0120
0.0125
```

The frozen Gate-7 value `0.012` lies in the middle of that interval. At those three points the matched/separated write ratio rises from about `7.77x` to `25.68x`, while the separated later-route change remains below `1.76%`.

This is useful but **not** a broad robustness result. The threshold window is finite and narrow. Below it, the separated event writes too much; above it, even the matched event becomes too weak to clear the required 5% route-change criterion.

## Write-scale sweep

Keep the threshold fixed at `0.012` and vary `eta` from `100` to `1600`, a 16-fold range. Absolute memory strength changes, as it should, but route selectivity survives:

```text
eta 100     matched/separated route change  12.43x
eta 800                                      9.67x
eta 1600                                     8.11x
```

Across the full sweep the matched/separated later-route selectivity stays between `8.11x` and `12.43x`.

So Gate 8 changes the status of Gate 7 from

```text
one flattering hand-picked point
```

to

```text
a real but parameter-sensitive operating region
```

The deterministic data are in [`results/gate8_attack_receipt.json`](results/gate8_attack_receipt.json). Run:

```bash
python gate8_attack.py
```

## What should be attacked next

The next test should stop relying on a threshold that makes A and B individually write zero. Allow all worlds to change and use exact inclusion/exclusion to isolate the persistent structural component attributable specifically to the nonlinear A+B encounter.

After that, move to the audited `Operaattori` morphology: local length/diameter edits, pure pose as a required null, operator tangents, and bounded readout. The point is no longer to add another decorative mechanism. It is to find out whether the causal write and global transfer story survives a real dendritic tree.
