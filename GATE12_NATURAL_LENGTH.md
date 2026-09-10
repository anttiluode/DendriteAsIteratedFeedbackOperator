# Gate 12 — Natural length is a grounded fixed point, not a wavelength slogan

The intuition was simple:

> If geometry compiles the dendritic operator, can growth stop at a length where the useful mode no longer improves under nearby length edits?

The first implementation found something more important than an immediate yes.

## 1. Length really is an operator variable

For one terminal dendritic compartment, Gate 12 changes a dimensionless length scale `s` while holding radius and channel densities fixed. In the toy this means

```text
g_axial -> g_axial / s
C        -> s C
g_leak   -> s g_leak
g_q      -> s g_q
```

So changing length does not act like multiplying one synaptic weight. It changes axial coupling, membrane load and the complex frequency response together.

At the 1% bounded-observer reference fixed point (`s = 0.62308`), the target transfer phase shifts by `0.14888 rad`, and the five-tone response still has a `2.55%` residual after fitting the best possible single complex scalar multiplier. Geometry is changing the operator, not just its gain.

## 2. The naive natural-length idea fails

If the growth rule maximizes only

```math
P_{\mathrm{pure}}(L)=
\frac{P_{\mathrm{target}}(L)}
     {P_{\mathrm{target}}(L)+P_{\mathrm{off}}(L)},
```

then the local search runs all the way to the imposed upper length bound.

```text
baseline target fraction       0.976354
scale 50 target fraction       0.984565
```

That looks like improvement until absolute transmission is inspected:

```text
target gain at scale 50 / baseline target gain = 0.002335
```

The branch became purer mostly by becoming almost disconnected.

This is the **silence trap**:

> a mode can dominate because everything else died.

That is directly analogous to the caution already learned in `SighImageSuper`: the final invariant mode can be perfectly dominant while the absolute signal is vanishingly small. Modal dominance by itself is not a useful answer.

## 3. Ground the geometry in a bounded observer

A bounded downstream observer has an absolute noise floor. Gate 12 therefore evaluates

```math
Q_\sigma(L)=
\frac{P_{\mathrm{target}}(L)}
     {P_{\mathrm{off}}(L)+\sigma^2},
```

where the noise floor is expressed as a counted fraction of the baseline target power.

Now infinite isolation is bad: once the target falls below the observer floor, the score collapses.

The local length update is deliberately tiny and explicit. For log-step `d`, evaluate

```text
L exp(-d), L, L exp(+d)
```

and move to whichever has the larger grounded score. Reduce `d` when the current length wins. A fixed point means that, at the current resolution, neither a slightly shorter nor slightly longer segment is better.

This is the structural analogue of the Sigh fixed-point idea:

```text
state fixed point:      x* = F(x*)
structure fixed point:  L* = G_sigma(L*)
```

No claim is made that a biological dendrite literally runs this optimizer.

## 4. Finite fixed points appear — but they are observer-relative

GitHub Actions produced:

| observer noise floor (fraction of baseline target power) | fixed-point length scale | target gain / baseline | target fraction |
|---:|---:|---:|---:|
| 0.1% | 1.5333 | 0.7291 | 0.97793 |
| 0.3% | 1.0717 | 0.9583 | 0.97663 |
| 1% | 0.6231 | 1.2423 | 0.97435 |
| 3% | 0.2284 | 1.5132 | 0.97092 |
| 10% | hits minimum tested length | 1.5321 | 0.97061 |

For the first four noise floors the broad sweep and the local update independently find interior optima. At 10% noise the optimum leaves the tested range toward a shorter segment.

So there is **not one universal natural dendritic length** in this toy.

There is a family

```math
L^*(\text{target},\text{morphology},\text{membrane},\text{observer floor}).
```

That is a better result than the original wavelength story. The structure is tuned not merely to an intrinsic eigenmode, but to a mode that must remain recoverable by a bounded downstream system.

## 5. What this says about growth

A plausible abstract growth loop is now:

```text
current geometry
      ->
frequency-dependent dendritic operator
      ->
mode mixture reaching AIS
      ->
bounded consequence/readout
      ->
local shorter/current/longer comparison
      ->
geometry update
      ->
repeat until local operator score stops changing
```

The stopping condition is not "phase error equals zero" and not "branch length equals one wavelength." In this gate it is simply

```math
G_\sigma(L^*)=L^*.
```

That is the clean meaning of **natural length operator** we can currently defend.

## 6. Why this matters for the larger line of repos

`SighImageSuper` supplied the warning that a surviving mode can be ungrounded.

`DendriteAsIteratedFeedbackOperator` supplied a physical geometry whose length changes the transfer operator.

Gate 11 supplied an event/receipt/write/replay loop.

Gate 12 now says that a useful structural write cannot optimize purity in isolation. It has to preserve a distinction strongly enough that a bounded observer can still recover it.

So the revised loop is

```text
PURIFY
  -> COMMIT
  -> RETURN RECEIPT
  -> PROPOSE STRUCTURAL CHANGE
  -> ASK: does the next bounded readout become more useful?
  -> KEEP / REVERSE
  -> repeat
```

This is much closer to growth as **operator search under consequence** than to ordinary scalar-weight learning.

## Claim boundary

Established only in this numerical toy:

- local segment length changes multiple cable terms together and changes complex transfer non-scalarly;
- purity-only geometry optimization has a pathological isolation solution;
- adding an explicit bounded-observer noise floor creates finite local length optima over a substantial tested range;
- those optima depend strongly on the observer floor.

Not established:

- that biological dendrites optimize this score;
- that a bAP supplies the sign of a length update;
- that dendritic branches grow toward literal standing-wave wavelengths;
- that the toy's dimensionless length scales map directly to micrometres;
- that the same result survives a real reconstructed morphology.

The next serious gate is therefore obvious: **apply the same length derivative / bounded-readout objective to the real Operaattori morphology, where length and diameter are actual geometric parameters rather than toy scale factors.**
