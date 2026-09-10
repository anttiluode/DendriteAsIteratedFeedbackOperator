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

# Gate 5 — active membrane state changes the spectrum

The passive cable was intentionally boring: it is a stable low-pass operator.

Gate 5 adds a minimal delayed restorative current on the distal half of branch B,

```math
C\dot V=-GV-g_rSw+I,
\qquad
\tau\dot w=S^TV-w.
```

This is **not** fitted as a quantitative HCN model. It is the smallest active-conductance motif needed to test the operator idea honestly.

The result changes qualitatively:

```text
passive peak / DC gain = 1.000
passive peak omega     = 0

active peak / DC gain  = 2.597
active peak omega      = 0.0791
complex poles          = 6
largest Re(pole)       = -0.0434
```

So the same cable topology can move from a purely low-pass transfer to a stable frequency-selective operator when active state is added.

That is the first place where the resonant-cavity picture becomes biologically useful without claiming the dendrite is literally an acoustic cavity:

> **morphology supplies the spatial constraints; active conductance supplies additional dynamical coordinates; together they determine which temporal modes are preferentially transmitted.**

---

# Gate 6 — the present dendritic state changes the operator

Gate 6 adds two local voltage-dependent NMDA-like conductances on the distal branch. The voltage gate is intentionally qualitative rather than a fitted receptor model.

Two coincident inputs produce a soma interaction residual whose norm is

```text
matched coincidence   1.08268e-4
15-time-unit separated 5.03242e-5
matched / separated   2.151
```

The more important test freezes the **same external synaptic gates** at `t = 8` and evaluates the instantaneous Jacobian twice:

```math
J(V,t)=\frac{\partial \dot V}{\partial V}.
```

Once at the actual A+B voltage state, and once at `V=0`.

Only the voltage-dependent internal state differs.

The Jacobian change has

```text
||Delta J|| = 0.010777
rank        = 2
```

because only two local nonlinear conductances are active.

Yet a tiny B-site probe sees a different soma transfer:

```text
zero-state linearization  0.145205
A+B-state linearization   0.171316
ratio                     1.17982
```

The complete A+B trajectory remains locally stable in this toy; the largest sampled real Jacobian eigenvalue is `-0.01299`.

This is the line I wanted the repo to reach:

```text
morphology
    -> baseline operator

current dendritic state
    -> different instantaneous operator

local nonlinear encounter
    -> local Delta J

same later perturbation
    -> different global response
```

So the branch is not merely carrying a state through a fixed filter.

\[
\boxed{\text{the state participates in defining the filter}}
\]

That is substantially closer to the object we kept circling in Sigh, Jello, the nonlinear-fluid machine, and the Geometric Neuron work.

The frozen active receipt is [`results/active_receipt.json`](results/active_receipt.json). Run it with:

```bash
python gate5_gate6.py
```

---

# Gate 7 — coincidence writes a persistent local constraint

Gate 6 still only changed the operator while the nonlinear state was present.

Gate 7 finally crosses the next boundary.

The A+B training event drives a local midpoint voltage. A deliberately simple local plasticity rule integrates only suprathreshold activity,

```math
\Delta g
=
\eta\int
\max(V_{\rm write}-\vartheta,0)^2\,dt.
```

The rule has no access to the labels "A" or "B". It only sees the local voltage consequence.

Five matched-budget worlds are compared:

```text
W0     no event
WA     A alone
WB     B alone
WAB    simultaneous A+B
Wsep   A and B separated by 15 time units
```

For the frozen parameters:

```text
WA write       0
WB write       0
WAB write      0.126078
Wsep write     0.0102823
WAB / Wsep     12.26x
```

The written quantity changes one local axial conductance between compartments 23 and 24.

Then the important part happens:

```text
erase every fast voltage state
erase every synaptic gate
start recall from exact zero
inject the same unit-charge A cue
```

Only the changed physical constraint survives.

The later B-site peak is

```text
baseline          0.0864009
matched memory    0.0935157
separated memory  0.0871369
```

or

```text
matched route change    +8.23%
separated control       +0.85%
```

The soma peak, meanwhile, changes by only `-1.19%`.

So the same structural memory is about **6.94x more visible locally at the branch than at the soma** when relative changes are compared.

That is unexpectedly clean because it joins three of the older threads in one toy:

```text
local nonlinear coincidence
        ->
persistent constraint edit
        ->
fast state erased
        ->
same later cue
        ->
different branch route
        ->
soma mostly misses it
```

This is no longer merely "the current state changes the filter."

\[
\boxed{
\text{experience changed the dendritic constraint,
and the changed constraint altered a later signal}
}
\]

The frozen receipt is [`results/slow_write_receipt.json`](results/slow_write_receipt.json). Run:

```bash
python gate7_slow_write.py
```

The important limitation is equally explicit: the plasticity threshold and scale are engineered in this toy. Gate 7 proves that the complete causal loop can be made local and persistent in the compartmental substrate; it does **not** prove that this exact write rule is biological or that such a memory will self-organize without design.

---

# Next gates

1. **Threshold attacker.** Sweep the Gate-7 write threshold/scale and report the region where matched coincidence remains selective rather than choosing one flattering point.
2. **Causal write subtraction.** Allow A and B individual worlds to write too, then retain only the inclusion/exclusion component attributable specifically to their nonlinear encounter.
3. **Real morphology.** Replace the Y toy with the audited `Operaattori` reconstruction; use local length/diameter edits and retain pure pose as a required null.
4. **Operator tangent test.** Compare measured global transfer edits with the low-dimensional prediction from local morphology/conductance tangents.
5. **Bounded query.** Search for a stimulation/readout that reveals a real structural edit the soma trace cannot identify reliably.

The conceptual chain is now complete in the toy:

```text
structure
  -> propagation operator
  -> active state changes instantaneous operator
  -> coincidence produces a local slow write
  -> fast state disappears
  -> changed structure changes later propagation
  -> bounded observation determines whether memory is visible
```

That is almost exactly the architecture we kept reaching for separately in `SighImageSuper`, `Operaattori`, Active Dendrite, Jello, and the nonlinear-fluid machine.

---

## Claim boundary

This repo does **not** establish that biological dendrites implement holographic memory, that action potentials are standing waves, that nodes of Ranvier are resonant phase filters, or that this architecture beats transformers.

It now establishes seven small numerical/algebraic facts in transparent compartmental toys:

- repeated passive cable evolution selects slow modes;
- continued forcing produces the same grounded fixed-point structure as the Sigh recursion;
- morphologically different dendritic branches have temporal kernels that cannot be reduced to one scalar multiplier;
- a local axial-conductance edit produces an exact global rank-one change in the continuous resolvent;
- a delayed restorative active current turns the same branch into a stable frequency-selective operator;
- two local voltage-dependent conductances make the instantaneous Jacobian state-dependent, so the same tiny probe can have a different global consequence;
- a local coincidence-triggered conductance write survives exact fast-state erasure and alters a later cue's branch route, while the soma is much less sensitive to that edit.

The remaining question is no longer whether the pieces can be connected.

It is whether the same mechanism survives **parameter attacks and real morphology** without being hand-held into existence.
