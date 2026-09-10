# Gate 11 — PURIFY -> COMMIT -> RETURN RECEIPT -> WRITE -> REPLAY

This gate is the first attempt in this repo to assemble the dendritic purifier into a biological-looking event loop instead of studying each mechanism separately.

It is intentionally a small operator toy, not a fitted pyramidal neuron.

```text
mixed temporal cue
      ↓
quasi-active dendritic purifier
      ↓
AIS boundary
      ↓ event
axon-side commitment
      +
post-event return through dendrite
      ↓
local receipt × local eligibility × write permission
      ↓
persistent local edge edit
      ↓
fast state erased
      ↓
same cue sees a different operator
```

## 1. SELECT / PURIFY

The existing `QuasiActiveDendrite` is retained.  One AIS compartment is appended to the soma through reciprocal axial conductance.  We find the non-zero frequency `omega_star` that maximizes distal-tip -> AIS transfer, then present five equal-amplitude temporal tones.

The target tone begins as only 20% of input energy.  Gate 11 records how strongly the dendritic operator enriches it at the AIS.

## 2. Small-signal reciprocity is the control

For a frozen physiological state,

```math
Y(\omega)=G+i\omega C+D_q(\omega)
```

is complex symmetric in this toy.  Therefore

```math
H(\omega)=Y(\omega)^{-1}
```

is also complex symmetric and

```math
H_{AIS,j}=H_{j,AIS}.
```

This is the baseline.  The gate explicitly checks `H == H.T` numerically.

## 3. COMMIT / RETURN: the receipt is history-conditioned, not an adjoint

An AIS event is represented by a transient post-spike conductance state: strong AIS refractory/AHP-like conductance plus a smaller somatic shunt.

The backward event therefore does **not** traverse the same frozen operator that carried the forward cue.

```text
forward small signal     H_pre
return after event       H_post
```

Define

```math
\delta r_j = [H_{post}]_{j,AIS}-[H_{pre}]_{j,AIS}.
```

Gate 11 calls `delta r` the **receipt**.

This is deliberately *not* called an adjoint and it is not proof of biological non-reciprocity.  Each frozen state is still reciprocal.  The asymmetry comes from history: the forward and return traversals occur in different channel states.

## 4. BIND / WRITE

The forward cue leaves a local eligibility magnitude `e_j`.  The event leaves a return-receipt magnitude `|delta r_j|`.

The first closed-loop rule is intentionally simple:

```math
score_j = e_j |\delta r_j|.
```

Only branch-A edges are candidates.  The edge with the largest local overlap receives

```math
\Delta g = \eta\, permission\, score_j.
```

`permission` is a separate consolidation gate.  It is the placeholder for the lesson from `JelloBrain`: a successful/stable route should not automatically reinforce itself forever.  Future gates can replace this bit with measured failure, surprise, relevance or delayed consequence.

The important distinction is

```text
receipt != learning
```

The receipt is a fast consequence signal.  Learning is the slower persistent parameter edit.

## 5. One local write -> one global structured operator change

The local axial edit has form

```math
Y' = Y + \Delta g\,bb^T.
```

So at fixed frequency the resolvent change is rank one by Sherman-Morrison.  Gate 11 records `s2/s1` of the measured `Delta H` as an audit.

This is the dendritic version of the same object that appeared in `Kompressori` and `resonant-graph-with-a-nonlinear-fluid-interior`:

```text
large response operator
      +
small local material change
      ↓
low-rank but distributed global deformation
```

## 6. REPLAY after fast-state erasure

The identical cue is then evaluated from a fresh frequency-domain state.  No previous voltage, gating trajectory or bAP state is carried into the replay.

Only the persistent edge edit remains.

The falsifiable question is:

> Does the same cue now reach the AIS differently because the preceding event changed the operator, rather than because fast voltage was left behind?

A `permission=0` control must produce exactly zero persistent edit and therefore no replay change.

## What this gate establishes if CI passes

Only the following narrow chain:

```text
mode-selective dendritic transfer
      ->
reciprocal small-signal forward/reverse control
      ->
AIS-state-conditioned return residual
      ->
local eligibility/receipt overlap
      ->
local persistent conductance edit
      ->
global low-rank transfer change
      ->
changed response to the same later cue after fast-state wipe
```

It does **not** establish that:

- biological bAPs compute an adjoint or gradient;
- the toy refractory conductances quantitatively model Nav/Kv/AHP kinetics;
- the eligibility rule is STDP;
- chandelier cells implement the permission bit;
- a real neuron performs this exact loop;
- the mechanism improves learning on a useful task.

Those are later attacks, not assumptions.

## Next attacks

The next useful gates are already forced by the architecture:

1. replace the event-state switch with a dynamical AIS spike model and measure the return waveform in time;
2. compare basket-like somatic shunt with chandelier-like distal-AIS shunt under matched conductance budget;
3. replace the permission bit with delayed failure/relevance and attack unconditional Hebbian reinforcement;
4. run many writes and measure when local low-rank edits add versus interfere;
5. transfer the loop to the real Operaattori reconstruction.
