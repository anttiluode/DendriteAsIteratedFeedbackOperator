# Gate 6 — nonlinear sieve, then the attacker

Gate 5 established the narrow result this repo now calls **eigenmode purification**: a deliberately quasi-active branch has a non-zero pass band, and an equal-energy five-tone input becomes strongly enriched in the branch-selected temporal mode at the soma.

Gate 6 asked the next question:

> Does that distributed mode selection determine which local activity enters a nonlinear regime?

## First mechanism test

Two sources on branch A were calibrated so that, at the branch's preferred frequency `omega = 0.125557`, their linear phasors contribute equally at a proximal hotspot.

The hotspot contains a bounded regenerative current with a knee at `|V| = 0.06`. It is intentionally generic and is **not** an NMDA kinetic model.

The matched pair produces a linear hotspot peak of `0.0800`; quadrature gives `0.05657`; either single gives `0.0400`.

After adding the nonlinear hotspot, the soma residual relative to the corresponding linear system is:

| condition | nonlinear residual RMS |
|---|---:|
| matched phase + branch carrier | `4.2039e-3` |
| quadrature | `2.8811e-3` |
| antiphase | ~`4e-20` |
| high off-band carrier | `6.93e-10` |
| spatially separated second input | `1.60e-6` |
| single input | `1.48e-6` |

So the matched condition beats the hardest control, quadrature, by `1.459x` and strongly beats the easier controls.

That demonstrates a **coherent coincidence -> local nonlinear event -> changed soma response** mechanism.

It does **not** yet demonstrate that dendritic resonance caused the selectivity.

## Gate 6b attacker 1 — move the nonlinear knee

The obvious criticism is that the knee was deliberately placed between the matched and quadrature linear amplitudes.

So nothing except the knee was changed. Carrier, phases, amplitudes, geometry and quasi-active substrate were frozen.

The matched/quadrature nonlinear-residual ratio was:

```text
knee 0.045 ->    1.151x
knee 0.050 ->    1.197x
knee 0.055 ->    1.278x
knee 0.060 ->    1.459x
knee 0.065 ->   41.822x
knee 0.070 ->  261.534x
knee 0.075 -> 1317.572x
```

So the result is not confined to the single chosen threshold. The advantage is modest at low knees and becomes enormous once quadrature falls below the regenerative regime.

This attacker therefore **survives**.

## Gate 6b attacker 2 — remove the resonant membrane

This is the important one.

All quasi-active conductances were set to zero. The raw protocol was kept **exactly identical**:

- same carrier `omega`;
- same source amplitudes;
- same relative phase;
- same source locations;
- same hotspot;
- same nonlinear current;
- same knee.

No passive model re-tuning was allowed.

The passive ablation still produced:

```text
matched linear hotspot peak       0.06633
matched nonlinear residual RMS    0.002927
quadrature residual RMS            1.93e-8
matched / quadrature               1.52e5
```

The active resonant substrate made the matched hotspot about `1.206x` larger and the matched nonlinear residual about `1.436x` larger, but **resonance was not necessary for the coincidence selectivity**.

In fact, under this particular raw protocol the passive matched/quadrature ratio is larger, because the phase calibration happens to cancel quadrature more completely in the passive substrate.

Therefore the strong version of Gate 6 is **not established**:

> We cannot say that the resonant dendritic purifier caused the nonlinear selection in this two-tone coincidence experiment.

That is not a reason to discard Gate 5. It tells us the next experiment has to force the purifier to do indispensable work.

## What Gate 6 actually establishes

The honest result is now:

```text
Gate 5:
    active distributed dendritic operator
    -> non-zero temporal mode purification          ESTABLISHED IN TOY

Gate 6:
    coherent same-branch activity
    -> local regenerative nonlinear event           ESTABLISHED IN TOY

Gate 6b:
    resonance is necessary for that nonlinear event NOT ESTABLISHED
```

The two mechanisms coexist in the same toy, but Gate 6 has not yet causally linked them.

## The required next gate

Do not give the nonlinearity a clean single carrier that a passive cable can also sum coherently.

Give both substrates the **same broadband/multimode input**, for example the five equal-energy components from Gate 5. The target mode begins as only 20% of input energy.

Then use exactly the same local hotspot nonlinearity and ask:

```text
ACTIVE QUASI-RESONANT BRANCH
broadband mixture
    -> target mode purified
    -> hotspot crossing
    -> nonlinear residual

PASSIVE ABLATION
same broadband mixture
    -> target remains weak / low-pass mixture
    -> ? hotspot crossing
    -> ? nonlinear residual
```

No re-tuning after ablation.

That is the experiment that can establish the stronger sentence:

> **The dendritic operator purified a mode, and that purified mode was what triggered the nonlinear event.**

Until that test is passed, the repo should keep the two claims separate.
