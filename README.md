# Dendrite As Iterated Feedback Operator

## Latest — Gate 11 closes the first biological-looking loop

The repo now goes beyond the isolated **eigenmode purifier** and connects it to an AIS-like event boundary, an event-conditioned return receipt, a local structural write, and replay after fast-state erasure.

```text
PURIFY
mixed temporal cue -> branch-selected non-zero mode

COMMIT
purified state reaches an AIS-like boundary event

RETURN
post-event state sends an AIS-originating receipt back through the same morphology

BIND / WRITE
local forward eligibility × local receipt × consolidation permission
        -> one persistent local conductance edit

REPLAY
erase fast state exactly
        -> same cue sees the changed operator
```

The first frozen GitHub Actions receipt is [`results/gate11_ais_receipt.json`](results/gate11_ais_receipt.json); the full mechanism and claim boundary are in [`GATE11_AIS_RECEIPT.md`](GATE11_AIS_RECEIPT.md).

First deterministic numbers:

```text
input target-mode energy                    20.00%
AIS target-mode fraction after purifier     97.64%
small-signal reciprocity error              1.76e-16
post-event return / small-return change     77.63%
local write Delta g                         1.576e-3
Delta H s2 / s1                             2.43e-12
same-cue AIS gain change after fast wipe    2.99e-4 relative
no-write-permission operator change         0
```

The important correction is explicit:

> **The receipt is not learning.**

The frozen pre-event and post-event operators are each reciprocal. The forward cue and return event differ because they occur in different channel states. Learning happens only when the local receipt meets local forward eligibility *and* a separate consolidation-permission signal, producing a slower persistent operator edit.

This is not a fitted AIS/bAP/STDP model and it does not claim biological backpropagation computes an adjoint or gradient. It is the first executable version here of:

> **PURIFY -> COMMIT -> RETURN CONSEQUENCE -> WRITE -> NEXT INPUT SEES A DIFFERENT OPERATOR.**

---

This repo started with a cautious bridge from `SighImageSuper` to passive cable theory. Gate 5 tests the stronger purifier picture directly:

> **Can a dendritic branch behave as an eigenmode purifier: taking a mixed temporal signal, rejecting most components, and delivering a branch-selected non-zero mode to a bounded soma readout?**

In the current numerical toy, **yes — operationally, in that narrow filtering sense.**

The claim is not that dendrites are acoustic cavities, RF waveguides, or hidden Helmholtz resonators. The useful object is the **distributed operator** compiled by morphology and membrane kinetics.

```text
local geometry + membrane kinetics
        ->
frequency-dependent propagation operator
        ->
mode-selective dendritic transfer
        ->
bounded soma / branch readout
```

The passive cable remains the control. The new active model is a phenomenological restorative quasi-active membrane, not a fitted HCN model.

---

# Gate 5 — non-zero mode purification

The active branch adds one restorative state `z` per dendritic compartment:

```math
C\dot V=-GV-g_q z+I,
```

```math
\tau_q\dot z=V-z.
```

For sinusoidal drive this contributes

```math
Y_q(\omega)=\frac{g_q}{1+i\omega\tau_q},
```

so the voltage transfer is

```math
H(\omega)=
\left[
G+i\omega C+
\mathrm{diag}\,Y_q(\omega)
\right]^{-1}.
```

This enlarged state system can have stable complex poles and a genuine non-zero pass band.

The deterministic GitHub Actions receipt is frozen in [`results/active_receipt.json`](results/active_receipt.json).

## Passive control: blur / diffusion

The matched passive morphology peaks at DC on both distal-to-soma channels:

```text
branch A peak omega = 0
branch B peak omega = 0
```

That is the old result. Repeated passive evolution selects the slowest decaying direction, but it is fair to call that **diffusive mode selection**, not the interesting purifier we were after.

## Active branches: two different non-zero pass bands

With quasi-active restorative membrane:

```text
branch A
    peak omega       0.125557
    peak / DC gain   29.91x
    Q-like           1.654

branch B
    peak omega       0.068077
    peak / DC gain   34.31x
    Q-like           1.388
```

The branch peak frequencies differ by a factor of

```text
1.844x
```

although the two branches terminate at the same soma.

So the branch is not merely a scalar coefficient attached to an input. It is a **tuned temporal transfer operator**.

## The actual purifier test

Give the distal tip five equal-energy temporal components. The target resonant component therefore starts with only

```text
20% of input energy.
```

At the soma, after propagation through the active dendritic operator:

```text
branch A target fraction: 97.61%
branch B target fraction: 96.03%
```

The same target frequencies through the passive control account for only

```text
branch A passive target fraction: 0.64%
branch B passive target fraction: 4.36%
```

So in this toy the dendritic channel really is doing the thing the original visual intuition suggested:

```text
mixed temporal input
        ->
branch operator
        ->
most components strongly rejected
        ->
selected non-zero mode dominates the bounded output
```

That is why the repo now uses **eigenmode purifier** as an operational hypothesis rather than as a metaphor.

Important: this is still linear resonant filtering. No information is created and no nonlinear winner-take-all process is implied.

## Complex poles are actually present

The full `V + z` state operator is stable:

```text
maximum real pole part = -0.02544
```

and contains

```text
12 complex poles.
```

The least-damped complex pole measured in the receipt is

```text
lambda = -0.03997 +/- i 0.05469
```

so this is categorically different from the passive RC control, whose modal generator has only real decays.

Adding a gating variable does **not** mean the dendrite has literally become an acoustic wave equation. It means the distributed electrical system now has oscillatory state-space modes and resonant transfer.

---

# Local structural edits still make global low-rank operator edits

The older passive gate gave an exact Sherman-Morrison result. For an axial edge edit,

```math
G\rightarrow G+\delta g\,bb^T,
```

the full resolvent changes by rank one.

Gate 5 shows the same structural idea survives active frequency selectivity.

At fixed `omega`, changing one local quasi-active conductance is one diagonal admittance edit. Numerically:

```text
Delta H effective rank        1.0000000000009
top singular fraction         0.999999999999946
s2 / s1                       1.05e-14
entries > 2% peak change      23.52%
```

So a **single local membrane change** can create a distributed global edit in the frequency-specific dendritic transfer operator.

This is the bridge back to `resonant-graph-with-a-nonlinear-fluid-interior` and `Kompressori`:

> **local physical change -> global structured operator change**

without storing or editing a dense matrix coefficient-by-coefficient.

---

# Why this is a stronger SighImageSuper connection

`SighImageSuper` began with

```math
x_{n+1}=Ax_n
```

and showed that the operator determines which differences disappear quickly and which persist.

The first dendrite gate copied only that passive logic:

```text
many spatial modes
    -> repeated cable decay
    -> slowest mode remains last
```

Gate 5 moves to the more interesting object:

```text
mixed temporal modes
    -> distributed resonant operator
    -> one non-zero branch-selected band survives the route
    -> bounded soma sees a purified answer
```

That is much closer to the Berglund visual intuition: geometry/material constraints determine **which global response is easy for the system to support**.

The crucial distinction is that we are not storing the selected waveform as a template.

```text
The structure makes that response natural.
```

Or in the language recurring across the repos:

> **Structure compiles the operator; the operator decides which modes can live.**

---

# Earlier gates

The current repo still preserves the passive controls because they tell us exactly what Gate 5 added.

## G0 — stable recurrent cable

A 25-compartment Y cable has passive step spectral radius

```text
0.9791865
```

so every autonomous passive mode decays.

## G1 — passive Sigh-style mode selection

A distal pulse begins with modal effective dimension

```text
10.005
```

and falls under silent iteration to

```text
step 20    3.741
step 80    2.830
step 240   1.514
step 600   1.019
```

with late absolute cosine `0.99863` to the slowest passive mode.

Useful control; not yet a resonant purifier.

## G2 — continued forcing preserves cue identity

The Sigh recursion

```math
V_{n+1}=\alpha q+(1-\alpha)PV_n
```

has fixed point

```math
V^*=\alpha[I-(1-\alpha)P]^{-1}q.
```

The numerical error is about `1.1e-15`. Two distal cues have nearly orthogonal grounded fixed points (`cos ~1.84e-6`) but converge toward the same passive late mode after forcing is removed (`cos ~0.99909`).

## G3 — branch != scalar weight

Equal charge injected into the two distal tips yields soma kernels peaking at different times:

```text
branch A: step 79
branch B: step 93
```

After fitting the best scalar multiplier, `18.23%` of one waveform remains unexplained.

A branch therefore carries delay/filter geometry that a scalar weight cannot represent.

## G4 — one axial edit -> rank-one global resolvent edit

For

```math
H(s)=[G+sC]^{-1}
```

and one edge edit `G -> G + delta_g bb^T`, Sherman-Morrison gives

```math
\Delta H
=-\frac{\delta g Hbb^TH}
        {1+\delta g b^THb}.
```

The numerical effective rank is exactly `1.0` to roundoff.

---

# What Gemini's stronger picture gets right — and what still needs testing

The useful criticism of the passive version was correct:

> **A passive RC cable only proves diffusive spectral decay. It does not demonstrate selective non-zero resonance.**

Gate 5 fixes that.

A few stronger statements should still remain hypotheses rather than conclusions:

- Active conductances can generate resonance and complex poles without making the dendrite literally a hyperbolic acoustic waveguide.
- NMDA gives voltage-dependent cooperative nonlinearity, but whether it sharpens *the particular resonant mode selected by the branch* must be measured rather than assumed.
- Branch points and spine necks create impedance filtering, but ordinary dendrite theory does not imply that every neck is a cavity-style tunneling resonator.
- The Gate 11 return receipt is history-conditioned because the event changes channel state; it is not an exact adjoint.
- A receipt is not enough for learning: persistent consolidation requires a slower write rule and a separate permission/consequence condition.

The right thing is therefore not to retreat from the eigenmode-purifier thought. It is to **make each stronger part executable and attack it with controls.**

---

# Next gates

The first closed loop now exists. The next attacks should make it less engineered:

1. replace the Gate 11 event-state switch with a dynamical AIS spike and measure the time-domain bAP waveform;
2. compare basket-like somatic shunt with chandelier-like distal-AIS shunt under equal conductance budget;
3. replace the consolidation-permission bit with delayed failure/relevance and attack unconditional reinforcement;
4. stream many events and measure when local low-rank edits remain additive versus interfere;
5. transfer the closed loop to the audited Operaattori reconstruction.

---

# Run

```bash
python -m pip install -e .[dev]
pytest -q
python experiment.py
python active_experiment.py
python gate11_ais_receipt.py
```

CI runs the full gate sequence on Python 3.10 and 3.12.

See [`THEORY.md`](THEORY.md), [`PURIFIER_CAUSALITY.md`](PURIFIER_CAUSALITY.md), and [`GATE11_AIS_RECEIPT.md`](GATE11_AIS_RECEIPT.md) for equations and claim boundaries.

---

## Claim boundary

Established by the numerical toys in this repo:

- passive cable evolution performs diffusive modal selection;
- continued forcing gives a cue-specific resolvent response;
- different branches have non-scalar temporal kernels;
- one local passive edge edit creates an exact rank-one global resolvent edit;
- a restorative quasi-active extension creates stable complex poles and non-zero distal-to-soma pass bands;
- from five equal-energy input tones, the chosen resonant component can dominate the soma/AIS output at about `96-98%` in deliberately tuned toys;
- one local active-conductance or axial-conductance edit can make a distributed rank-one transfer edit at fixed frequency;
- Gate 11 verifies reciprocal small-signal forward/reverse transfer, a distinct event-conditioned return through a changed fast conductance state, local receipt/eligibility-gated persistent writing, and changed replay after fast-state erasure.

Not established:

- that real dendrites routinely achieve this degree of purification;
- that HCN alone implements the exact toy parameters;
- that biological dendrites are acoustic/RF waveguides;
- that NMDA automatically performs mode competition;
- that biological bAPs compute an adjoint, error gradient, or exact credit assignment signal;
- that chandelier cells implement the consolidation gate;
- that the mechanism gives useful learning, hardware efficiency, transformer equivalence, or brain-level computation.

Those are now experiments rather than metaphors.
