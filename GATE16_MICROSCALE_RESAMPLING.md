# Gate 16 — Microscale equal-budget response survives 2x morphology refinement

Gate 15 removed the apparent 20–1000 um support threshold once every support received the same total structural edit budget. One loophole remained: Gate 13's failed coordinate was only ~2 um, while Gate 15 started near 20 um.

Gate 16 bridges that decade directly and attacks numerical discretization at the same time.

## Protocol

The support selector is intentionally removed. Gate 15 showed that `|forward| x |return|` has not earned structural credit-assignment status. Gate 16 fixes the editable region at the terminal end of the driven soma-to-tip route and asks only about **support scale**.

```text
physical support widths        2, 4, 8, 16, 32, 64, 128 um
fixed total |Delta L| budget   0.2 um
observer floor                 1% of baseline target power
carrier                        175.868675 rad/s
morphology                     pinned human L2/3 cell 1125
```

Then every physical morphology edge is subdivided in two. The refined tree contains `25,027` dendritic nodes instead of `12,514`, but the physical endpoints and total cable are preserved.

```text
original path length   1242.169057660336 um
refined path length    1242.1690576603364 um
relative error         3.66e-16
```

This is a direct test of whether the small-support response is tied to arbitrary morphology-file sample intervals.

The compact frozen CI receipt is [`results/gate16_microscale_resampling_summary.json`](results/gate16_microscale_resampling_summary.json).

## Result — there is still no support turn-on

Original discretization:

```text
requested W     actual W      support nodes      Gamma_B / Q0
2 um             2.505 um          2             +0.02906%
4 um             3.229 um          3             +0.02904%
8 um             8.461 um          5             +0.02912%
16 um           17.150 um         10             +0.02912%
32 um           32.365 um         21             +0.02909%
64 um           63.708 um         40             +0.02908%
128 um         128.005 um         80             +0.02907%
```

The response is essentially flat from a two-edge ~2.5 um support to a ~128 um support when the total structural budget is held fixed.

So the remaining version of the proposed coherence threshold loses too:

> **In this model, useful local operator sensitivity extends to the smallest physical support resolved by the original morphology.**

There is no evidence here for a minimum coherent rewrite width between ~2 and ~128 um.

## The discretization attack

After subdividing every morphology edge in two, the same physical-width experiment gives:

```text
requested W     refined actual W    refined nodes    Gamma_B / Q0
2 um              1.512 um              3            +0.02927%
4 um              3.229 um              6            +0.02927%
8 um              8.461 um             10            +0.02932%
16 um            15.781 um             19            +0.02931%
32 um            32.365 um             42            +0.02930%
64 um            63.708 um             80            +0.02929%
128 um          128.005 um            160            +0.02929%
```

Curve agreement:

```text
cosine(original, refined)      0.9999998861
relative L2 difference         0.7256%
```

At matched widths the refined/original response ratio stays around `1.0066–1.0081`.

So the small-support effect is **not** explained by the original morphology discretization. A physical support represented by twice as many compartments gives essentially the same response curve.

The baseline transfer itself also changes only slightly under refinement:

```text
target gain refined/original             0.993302
five-tone best-scalar residual           0.000861
```

The refined model is therefore close enough to the original for this to be a meaningful numerical invariance check, while not being artificially identical.

## Every local winning move is still shorter

For all seven widths, in both discretizations, the best equal-budget perturbation is `shorter`.

This reinforces Gate 15's warning: the local improvement is behaving like a broad geometry sensitivity, not like a special resonance-matched support discovery mechanism.

## What this changes about Gates 13–14

The corrected reading is now:

```text
Gate 13
    one ~2 um coordinate had no finite natural-length fixed point under a broad scale search

Gate 14
    large coherent supports could have finite grounded optima when every segment
    received the same fractional scale change

Gate 15
    equalizing total structural budget removed the apparent 20–1000 um turn-on

Gate 16
    sensitivity continues essentially flat down to ~2 um and survives 2x numerical refinement
```

So Gate 13 should **not** be read as `tiny geometry edits do nothing`.

It showed something different: one microscopic coordinate, when allowed to wander over a large scale range, did not possess the kind of finite fixed point that the larger Gate-14 coordinates did.

That distinction matters:

```text
local sensitivity != existence of a finite structural fixed point
```

A tiny region can influence the operator while still lacking a self-contained optimum under the chosen growth rule.

## What survives of the "natural length" idea

The strongest surviving statement is no longer about a minimum editable support.

It is about **fixed points of coupled structural coordinates**.

Some extended coordinates can possess finite observer-grounded optima even though infinitesimal/local sensitivity exists everywhere. The interesting object may therefore be:

> **not the minimum scale at which geometry can change the operator, but the scale and organization at which many local sensitivities combine into a stable structural fixed point.**

That is a better question.

It shifts the next experiment from `where does sensitivity turn on?` to:

```text
when do many editable geometric degrees of freedom acquire a collective interior optimum?
```

That could still connect naturally to self-carving geometry, but it is a different claim from a Nyquist/coherence threshold.

## Claim boundary

The morphology is real, but the quasi-active membrane, observer objective, and geometry perturbation rule are phenomenological. Twofold edge refinement rejects one obvious discretization-artifact explanation; it is not a full numerical-convergence study. Nothing here establishes a biological minimum growth unit, a dendritic wavelength, or a structural credit-assignment mechanism.

The result is useful precisely because it narrows the idea:

> **The operator is sensitive to geometry all the way down. The unresolved question is how those local sensitivities organize into a stable grown form.**
