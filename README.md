# Dendrite As Iterated Feedback Operator

This repo started with a cautious bridge from `SighImageSuper` to passive cable theory. Gate 5 now tests the stronger picture directly:

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

# Latest result — Gate 5: non-zero mode purification

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

Important: this is still linear resonant filtering. No information is created and no nonlinear winner-take-all process is implied. The next gate asks whether voltage-dependent coincidence can turn this from a fixed resonant sieve into a **state-dependent** one.

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

The right thing is therefore not to retreat from the eigenmode-purifier thought. It is to **make each stronger part executable and attack it with controls.**

---

# Next gates

## Gate 6 — nonlinear resonant sieve

Put a local voltage-dependent coincidence element on the already-resonant branch.

Compare equal total input energy under:

```text
matched carrier + matched timing
matched carrier + wrong timing
wrong carrier + matched timing
spatially separated inputs
single-input controls
```

Measure whether the nonlinear residual is selectively largest when activity occupies the branch's own resonant pass band.

The important question is not merely "does NMDA make a bigger voltage?" It is:

> **Does local nonlinearity preferentially amplify the mode the distributed branch operator already selected?**

## Gate 7 — state-dependent Jacobian

Linearize before and during the nonlinear event:

```math
J(x,\theta)=\partial F/\partial x.
```

Then ask whether the selected mode changes the operator that selects the next mode.

## Gate 8 — slow write

Let successful coincidence alter one local conductance, erase fast state, and replay the identical broadband cue.

The target loop becomes:

```text
mode selected
    -> local nonlinear encounter
    -> local material edit
    -> global frequency-dependent resolvent changes
    -> future mode selection changes
```

That would be the dendritic analogue of the direct-fluid machine:

```text
wave -> encounter -> changed medium -> next wave travels differently
```

## Gate 9 — real Operaattori morphology

Transfer the whole experiment to the audited reconstructed morphology. The decisive test is whether real local length/diameter/conductance changes produce low-dimensional but global changes in **frequency-dependent** transfer and whether different real branches purify different temporal modes.

---

# Run

```bash
python -m pip install -e .[dev]
pytest -q
python experiment.py
python active_experiment.py
```

CI runs all of these on Python 3.10 and 3.12.

See [`THEORY.md`](THEORY.md) for the equations and claim boundaries.

---

## Claim boundary

Established by this numerical toy:

- passive cable evolution performs diffusive modal selection;
- continued forcing gives a cue-specific resolvent response;
- different branches have non-scalar temporal kernels;
- one local passive edge edit creates an exact rank-one global resolvent edit;
- a restorative quasi-active extension creates stable complex poles and non-zero distal-to-soma pass bands;
- from five equal-energy input tones, the chosen resonant component can dominate the soma output at about `96-98%` in this deliberately tuned toy;
- one local active-conductance edit still makes a global rank-one transfer edit at fixed frequency.

Not established:

- that real dendrites routinely achieve this degree of purification;
- that HCN alone implements the exact toy parameters;
- that biological dendrites are acoustic/RF waveguides;
- that NMDA automatically performs mode competition;
- that the mechanism gives useful learning, hardware efficiency, transformer equivalence, or brain-level computation.

Those are now experiments rather than metaphors.
