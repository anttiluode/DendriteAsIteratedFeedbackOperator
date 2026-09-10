# Gates 6–9 — from fixed purifier to state-dependent and persistent operator

The repository now contains **two parallel Gate-5 active constructions** that were developed at the same time and later reconciled rather than overwriting one another.

- `active_experiment.py` / `active_resonance.py` is the stronger fixed **eigenmode-purifier** result centered in the README: two branches acquire distinct non-zero pass bands and purify a five-component temporal mixture at the soma.
- `gate5_gate6.py` / `active_dendrite.py` is a smaller independent active-control branch whose main purpose is to lead directly into a state-dependent Jacobian and nonlinear coincidence test.

Their receipts are deliberately kept separate:

```text
results/active_receipt.json          fixed active purifier
results/gate5_gate6_receipt.json     small active control + Gate 6
```

The useful convergence is that both replace the passive picture

```text
fixed morphology -> low-pass decay
```

with

```text
morphology + membrane dynamics -> frequency-selective operator.
```

From there, Gates 6–9 ask the harder question: **can the operator change because of what is currently happening, and can that change outlive the fast state?**

## Gate 6 — state participates in defining the operator

Two local qualitative NMDA-like voltage-dependent conductances are added to the distal branch.

Matched A+B activity has a soma nonlinear interaction norm about `2.151x` a temporally separated control. More importantly, at the same external synaptic-gate state, the instantaneous Jacobian

```math
J(V,t)=\partial \dot V/\partial V
```

is evaluated once at the actual A+B voltage state and once at `V=0`.

The difference is local rank 2, but a small B-site probe sees a different global transfer:

```text
zero-state B -> soma transfer   0.145205
A+B-state B -> soma transfer    0.171316
ratio                           1.17982
```

So the fast state is not merely a signal moving through a fixed filter. In this nonlinear toy, **the fast state changes the filter through which the next perturbation travels.**

## Gate 7 — persistent local write

A local suprathreshold activity rule changes one axial conductance. After training, every fast voltage and synaptic gate is erased exactly and recall starts from zero.

The later identical A cue gives:

```text
B-site route change, matched memory      +8.23%
B-site route change, separated control   +0.85%
soma change, matched memory              -1.19%
```

The same structural memory is therefore about `6.94x` more visible locally at the branch than at the soma on this relative metric.

This completes the toy causal loop:

```text
nonlinear encounter
    -> persistent local constraint edit
    -> fast state erased
    -> same later cue
    -> different route.
```

## Gate 8 — parameter attacker

The Gate-7 write threshold is swept across 18 values. A strict selectivity criterion survives at three adjacent sampled values only:

```text
0.0115, 0.0120, 0.0125
```

That is a **finite but narrow** operating region, not broad robustness.

At fixed threshold, a 16-fold write-scale sweep keeps matched/separated later-route selectivity between `8.11x` and `12.43x`.

See [`GATE8.md`](GATE8.md).

## Gate 9 — remove the threshold and subtract individual writes

The privileged threshold is removed entirely:

```math
\Delta g=\eta\int V_{write}(t)^2dt.
```

Now A alone and B alone both write. Counterfactual inclusion/exclusion isolates the nonlinear persistent component:

```text
matched collision-specific write      0.09051
separated collision-specific write    0.05664
ratio                                  1.598x
```

Each world keeps its own total physical write. After exact fast-state erasure and the same later A probe, inclusion/exclusion is applied only to the measured responses:

```text
matched B-site chi norm       0.019409
separated B-site chi norm     0.012345
ratio                          1.572x
```

Across `eta = 10..100`, that causal-response selectivity remains `1.572x..1.618x`.

This is deliberately a **modest ~1.6x result**. It is more important than the larger thresholded Gate-7 number because individual writes are now explicitly allowed and accounted for.

See [`GATE9.md`](GATE9.md).

## Current claim boundary

The active conductances and plasticity laws are transparent qualitative toys, not fitted biophysical models. Inclusion/exclusion is an analysis across parallel counterfactual worlds, not an operation available to one neuron.

What the combined gates establish is narrower:

> A compartmental dendritic substrate can be a fixed frequency-selective operator; local nonlinear state can change its instantaneous Jacobian; a local activity consequence can change a persistent physical constraint; after fast-state erasure that changed constraint can alter a later route; and bounded observation can strongly affect whether that structural memory is visible.

The next serious test is **real morphology**, not another decorative toy gate.
