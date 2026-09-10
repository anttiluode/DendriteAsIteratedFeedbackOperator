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
