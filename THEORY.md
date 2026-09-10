# Theory note — from Sigh to cable theory

## 1. Same operator, different physical interpretation

`SighImageSuper` used repeated application of an image operator,

```math
x_{n+1}=Ax_n.
```

A passive dendrite after compartmentalization is

```math
C\dot V=-GV.
```

Over a finite interval,

```math
V_{n+1}=\exp(-\Delta t C^{-1}G)V_n.
```

So the common object is not "image" or "wave." It is repeated state evolution under a structured operator.

## 2. Why slow modes emerge

Because `G` is symmetric positive definite in the passive toy and `C` is diagonal positive, define

```math
B=C^{-1/2}GC^{-1/2}.
```

`B` is symmetric positive definite and has eigenpairs `(lambda_k, q_k)`. In energy coordinates `y=C^{1/2}V`,

```math
y(t)=\sum_k a_k e^{-\lambda_k t}q_k.
```

Modes with larger `lambda_k` disappear faster. Normalized late state tends toward the slowest mode with nonzero initial projection.

This is exactly the spectral-selection logic of the original recursive image loop, but with a cable generator rather than an FFT filter.

That passive result is a **control**, not yet the interesting "eigenmode purifier" claim. It is diffusive forgetting: all autonomous poles are real and the eventual normalized state is the slowest decaying direction.

## 3. Grounded recursion is a resolvent

For

```math
V_{n+1}=\alpha q+(1-\alpha)PV_n,
```

stable `P` gives

```math
V^*=\alpha[I-(1-\alpha)P]^{-1}q.
```

The persistent answer is therefore a resolvent response of the dendritic propagation operator.

That is close in mathematical form to the resonant cavity Green's function

```math
H(\omega)=[K-\omega^2I+i\gamma\omega I]^{-1},
```

although the physical spectra are different. Passive cable modes are decaying; resonant cavities have oscillatory modes.

## 4. A dendritic "weight" is a transfer function

With input matrix `B_in` and readout `C_out`,

```math
y(s)=C_{out}[G+sC]^{-1}B_{in}u(s).
```

A synapse therefore does not have only a scalar gain. Its effective influence includes location-dependent filtering, delay, interaction with morphology, membrane state, and active conductances.

A scalar weight is at best a low-dimensional summary of this transfer.

## 5. Local edit, global response

For an axial edge `(i,j)`, let

```math
b=e_i-e_j.
```

Changing its conductance by `delta_g` produces

```math
G'=G+\delta_g bb^T.
```

Let

```math
H=(G+sC)^{-1}.
```

Then

```math
H'
=
H-
\frac{\delta_g Hbb^TH}
     {1+\delta_g b^THb}.
```

Therefore

```math
\Delta H=H'-H
```

has rank at most one.

This is exact for the linear passive resolvent and is the same structural algebra as a one-edge edit in the resonant graph.

For multiple independent local edits,

```math
G'=G+U D U^T,
```

Woodbury implies a global transfer update whose rank is bounded by the number of edited local directions, before degeneracies.

That gives a concrete interpretation of low-rank global learning from sparse local physical change.

## 6. What active dendrites add

With voltage-dependent conductances and gating state `z`,

```math
\dot x=F(x,\theta,u).
```

The local effective operator is the Jacobian

```math
J(x,\theta)=\frac{\partial F}{\partial x}.
```

Now the operator is state-dependent. Two identical anatomical structures can route differently because their local voltage/gating states differ.

NMDA, sodium/calcium spikes, shunting inhibition, HCN and other conductances therefore fit naturally as mechanisms that move the dendrite from

```text
fixed structured operator
```

to

```text
state-dependent structured operator.
```

A slow plastic change then changes `theta` itself, producing a persistent operator edit.

That is the point at which the direct-fluid loop

```text
wave -> nonlinear encounter -> changed material -> future wave travels differently
```

and the dendritic loop become mathematically comparable without asserting they are physically the same substrate.

## 7. Gate 5: the stronger "eigenmode purifier" object

The first repo version stopped at passive spectral decay. Gate 5 deliberately asks the stronger question.

Attach one restorative quasi-active state `z` to each dendritic compartment:

```math
C\dot V=-GV-g_q z+I,
```

```math
\tau_q \dot z=V-z.
```

For sinusoidal drive, eliminating `z` gives the frequency-dependent local admittance

```math
Y_q(\omega)=\frac{g_q}{1+i\omega\tau_q},
```

and therefore

```math
H(\omega)=
\left[
G+i\omega C+
\operatorname{diag}Y_q(\omega)
\right]^{-1}.
```

This is intentionally phenomenological. It is a quasi-active restorative current, not a fitted HCN channel and not a claim that the cable equation has literally become an acoustic wave equation.

But mathematically it crosses the line the passive control could not cross: the enlarged state operator acquires stable complex poles and the distal-to-soma transfer can peak at a **non-zero** temporal frequency.

The GitHub CI receipt gives:

```text
passive branch A peak: omega = 0
passive branch B peak: omega = 0

active branch A peak: omega = 0.12556   peak/DC = 29.91x
active branch B peak: omega = 0.06808   peak/DC = 34.31x
peak-frequency ratio A/B = 1.844x
complex poles = 12
max real pole part = -0.02544  (stable)
```

So the two dendritic branches are no longer merely two differently blurred low-pass kernels. They are two differently tuned temporal operators.

The operational purification test is deliberately simple. Feed five equal-energy temporal components into a distal tip. The target resonant component begins with only `20%` of input energy. At the soma:

```text
branch A target fraction -> 97.61%
branch B target fraction -> 96.03%
```

For the matched passive cable the same target frequencies account for only `0.64%` and `4.36%` respectively.

That justifies a precise use of the word **purifier** in this toy:

> **A distributed dendritic operator can enrich one non-zero temporal mode from a broadband/equal-tone input while rejecting other components before the bounded soma readout.**

Nothing nonlinear is required for that first meaning; ordinary resonant filtering is enough. Nonlinearity becomes essential for the stronger next claim: state-dependent competition, coincidence-triggered gain, switching, and persistent operator rewriting.

There is also a useful structural continuation of Gate 4. At fixed `omega`, changing one local quasi-active conductance changes one diagonal admittance direction. The global transfer update remains numerically rank one:

```text
Delta H effective rank = 1.0000000000009
s2 / s1 = 1.05e-14
23.52% of global entries exceed 2% of peak change
```

So active resonance does not destroy the old local-edit/global-operator algebra. It makes it frequency-selective.

The stronger architecture is therefore now:

```text
local morphology + membrane kinetics
        ->
frequency-dependent global resolvent
        ->
non-zero mode selection / purification
        ->
local nonlinear coincidence
        ->
state-dependent Jacobian
        ->
slow structural edit
        ->
new future resolvent
```

That is substantially closer to the original Sigh/Berglund intuition than passive electrotonic decay alone.
