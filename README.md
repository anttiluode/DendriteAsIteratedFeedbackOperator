# Dendrite As Iterated Feedback Operator

A small falsification repo for a picture that became visible after `SighImageSuper` and `resonant-graph-with-a-nonlinear-fluid-interior` converged on the same object from opposite directions.

The claim is **not** that a dendrite is literally an acoustic cavity or that neuronal signaling is secretly standing-wave RF physics.

The narrower claim is stronger because it is ordinary mathematics:

> **A compartmental dendrite is already an iterated distributed operator. Morphology and conductance define that operator; repeated evolution selects a hierarchy of modes; synaptic location determines how those modes are entered; soma/branch readout determines which survive as observable consequences.**

That makes a dendritic branch much richer than a scalar weight.

---

## The bridge

For a passive compartmental cable,

```math
C\dot V=-G(\theta)V+Bu(t),
```

where

- `C` is compartment capacitance,
- `G(theta)` contains leak and axial conductances,
- `theta` is morphology/material structure,
- `B` says where current enters.

Over one time step,

```math
V_{n+1}=P_\theta V_n,
\qquad
P_\theta=\exp[-\Delta t\,C^{-1}G(\theta)].
```

That is the Sigh object.

Under silence,

```math
V_n=P_\theta^nV_0.
```

The operator gives different spatial directions different forgetting times. Repeated application progressively removes the faster directions and leaves the slowest supported modes.

So the analogy is not

```text
dendrite = waveguide because it looks like a channel
```

but

```text
structure -> operator -> repeated propagation -> modal persistence -> bounded readout
```

which is exactly the line that survived the audits in `SighImageSuper`.

---

# First deterministic receipt

Run:

```bash
python -m pip install -e .[dev]
pytest -q
python experiment.py
```

The frozen result is in [`results/receipt.json`](results/receipt.json).

## G0 — the cable really is a stable recurrent operator

The 25-compartment Y cable has exact passive step spectral radius

```text
0.9791865
```

so every autonomous mode decays.

This is important: the first bridge does **not** require resonance or a hidden oscillator. Ordinary passive cable dynamics already give the iterated-operator structure.

## G1 — Sigh-style mode selection happens in the cable

A unit distal pulse initially occupies many cable modes.

```text
modal effective dimension: 10.005
```

After repeated silent propagation:

```text
step 20    3.741
step 80    2.829
step 240   1.514
step 600   1.019
```

At step 600 the normalized state has absolute cosine

```text
0.99863
```

with the slowest eigenmode.

The amplitude is by then tiny. That is exactly the useful lesson from Sigh: **modal selection and memory amplitude are separate questions.** A dying system can become geometrically simple before it vanishes.

## G2 — grounding prevents everything becoming the same dying mode

Use the original Sigh recursion on the dendritic step operator:

```math
V_{n+1}=\alpha q+(1-\alpha)P_\theta V_n.
```

Its analytic fixed point is

```math
V^*=\alpha[I-(1-\alpha)P_\theta]^{-1}q.
```

The numerical iteration matches that expression to relative error

```text
9.99e-16
```

for `alpha = 0.08`.

Two different distal branch cues have almost orthogonal grounded fixed points:

```text
cosine = 1.84e-6
```

but if the fresh cue is removed and both states are allowed to decay for long enough, their normalized late states converge toward the same slow mode:

```text
late ungrounded cosine = 0.99909
```

So the clean neuronal translation of Sigh's `alpha` is **continued input / boundary forcing / reference drive**, not a magical biological "grounding parameter."

## G3 — a branch is not a scalar synaptic weight

Inject the same unit charge at the tips of the two morphologically different branches and observe only the soma.

The two impulse-response kernels peak at different times:

```text
branch A: step 79
branch B: step 93
```

Even after fitting the single best scalar multiplier from one kernel to the other, the relative waveform residual is

```text
18.23%
```

A scalar weight can match amplitude. It cannot in general match the branch's **temporal transfer kernel**.

This is the first precise form of the visual intuition:

> **the dendritic channel is itself part of the computation.**

## G4 — one local structural edit can make a global rank-one operator edit

At Laplace frequency `s`, the cable transfer is

```math
H(s)=[G+sC]^{-1}.
```

Changing one axial edge conductance by `delta_g` changes

```math
G\rightarrow G+\delta_g bb^T,
```

where `b=e_i-e_j`.

Sherman-Morrison therefore gives an exact rank-one change in the global resolvent:

```math
\Delta H
=
-\frac{\delta_g\,Hbb^TH}
       {1+\delta_g\,b^THb}.
```

The numerical receipt returns

```text
effective rank of Delta H = 1.000000
```

with the first singular value carrying effectively all singular mass.

This is the direct dendritic version of the result in `resonant-graph-with-a-nonlinear-fluid-interior`:

> **a local physical constraint edit can induce a global but low-rank change in how the whole object responds.**

That is a much cleaner candidate for a biological "weight update" than pretending every synapse is an isolated multiplier.

---

# What changed in the picture

The useful ontology is now:

```text
weight / parameter    = local physical constraint
branch operator       = cable propagation implied by all constraints
state                 = voltages / conductance states currently evolving
memory time           = modal / transient persistence of that operator
synaptic address      = where and when the state is injected
readout                = what soma / branch / downstream cell can observe
learning              = a lasting constraint change that alters later transfer
```

That is close to the resonant-fluid machine, but importantly it does **not** require the dendrite to be an acoustic resonator.

The shared object is the operator.

---

# Where the wave picture is biologically legitimate — and where it is not

Passive dendrites are well described at first order by cable equations and are predominantly low-pass/diffusive, not lossless waveguides. Active conductances can add frequency selectivity and genuine subthreshold resonance; dendritic HCN and related currents are established examples. Regenerative dendritic spikes add strong local nonlinearity.

So the next biologically serious extension is not to declare every dendrite a Helmholtz cavity. It is to add known active conductances and ask whether the same operator picture becomes richer:

```text
passive structure
    -> fixed P(theta)

active membrane state
    -> P(theta, x, gating)

local coincidence / NMDA / dendritic spike
    -> transiently different effective operator

slow plastic change
    -> persistently different operator
```

That is where the fluid result becomes relevant again.

---

# Important correction: nodes of Ranvier

The tempting "periodic resonant cavity / phase-matched repeater" description is too strong.

Nodes of Ranvier are active excitable gaps between myelinated passive cable segments. Myelin changes capacitance and resistance; concentrated voltage-gated channels regenerate the action potential. Internodal geometry strongly affects conduction velocity and safety factor, but ordinary saltatory conduction does **not** require the internode spacing to be tuned to a standing-wave wavelength for constructive interference.

So Ranvier is useful here as an example of **alternating passive propagation and active regeneration**, not as evidence that axons are RF resonator chains.

Keeping that distinction makes the dendrite idea harder to dismiss.

---

# Why this may matter for AI

A standard neural layer stores a large number of independent coefficients and then explicitly multiplies by them.

A physical or simulated operator substrate can instead store **constraints** whose collective propagation creates the effective dense map.

The resonant-graph repo showed this for cavity constraints. This repo shows the same algebra appears immediately in an ordinary compartmental cable.

The possible AI primitive is therefore not "replace matrix multiplication with waves." It is:

```text
local editable constraints
        ->
shared recurrent operator
        ->
input-dependent trajectories
        ->
small bounded readouts
```

with local edits producing structured global operator changes.

Whether that is computationally cheaper or more useful than conventional networks is completely unestablished. But it is now a testable architecture rather than an analogy.

---

# Next gates

1. **Active resonance.** Add an HCN-like gating variable and show how the branch transfer spectrum changes without pretending the passive cable already oscillates.
2. **NMDA / local nonlinear encounter.** Two equal inputs, matched versus separated in space/time; measure the non-additive soma and branch response.
3. **Jacobian operator.** Linearize the active dendrite around different local states and directly measure `J(theta, x)` as the instantaneous operator.
4. **Slow write.** Let a local coincidence change one conductance, erase fast state, then ask whether a later identical cue follows a measurably different route.
5. **Real morphology.** Replace the toy Y with the audited `Operaattori` morphology and test whether its local length/diameter tangents produce the predicted low-dimensional global transfer edits.

The last one is the strongest bridge: `Operaattori` already compiles morphology into a cable operator. `SighImageSuper` tells us how repeated operators organize persistence. The resonant-fluid machine tells us how local material change can alter later propagation. This repo is where those three statements become one equation.

---

## Claim boundary

This repo does **not** establish that biological dendrites implement holographic memory, that action potentials are standing waves, that nodes of Ranvier are resonant phase filters, or that this architecture beats transformers.

It establishes four small numerical/algebraic facts in a transparent compartmental model:

- repeated cable evolution selects slow modes;
- continued forcing produces the same grounded fixed-point structure as the Sigh recursion;
- morphologically different dendritic branches have temporal kernels that cannot be reduced to one scalar multiplier;
- a local axial-conductance edit produces an exact global rank-one change in the continuous resolvent.

Those are enough to justify the next experiment.
