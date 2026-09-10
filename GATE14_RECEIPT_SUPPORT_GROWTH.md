# Gate 14 — a natural length reappears when the structural coordinate has spatial support

Gate 13 failed cleanly on the real human morphology: the forward/return maximum named a single approximately `2 um` morphology interval, and changing that one interval barely changed the global operator. No tested bounded-observer objective had an interior length fixed point.

Gate 14 asks whether the mistake was not **length**, but **what we called one length variable**.

The Berglund-like object is a coherent channel. A MorphIO sample interval is merely a numerical point spacing. So this gate promotes the Gate-11-style return quantity from a point to a spatial field and coherently changes a connected region of the real soma-to-input cable.

The result is a scale transition.

## Fixed setup

Nothing in the membrane model was retuned after Gate 13.

The same pinned human L2/3 cell-1125 morphology is used:

```text
12,514 dendritic nodes including soma root
20.742 mm total dendritic cable
selected soma-to-tip route: 1,242.17 um
```

The same baseline quasi-active operator has

```text
omega* = 175.869 rad/s
baseline five-tone target fraction = 99.3765%
```

The baseline forward/return field is

```math
R_j = |H_{j,input}(\omega^*)|\,|H_{j,AIS}(\omega^*)|.
```

The peak is at the terminal input tip. That is an important warning: **this gate does not establish that forward × reciprocal return is a biologically meaningful credit locator.** It is only the fixed spatial field used to define the structural supports being tested.

Before any length sweeps, three support thresholds were fixed:

```text
R >= 0.9 Rmax
R >= 0.5 Rmax
R >= 0.1 Rmax
```

For each threshold, only the connected component containing the peak on the actual input path is changed coherently.

## Three structural scales

The resulting supports are very different physical objects:

```text
threshold   cable support    path fraction    receipt-mass fraction
0.9          42.31 um          3.41%               8.68%
0.5         356.60 um         28.71%              53.57%
0.1        1098.79 um         88.46%              98.10%
```

This immediately gives a stronger version of the Gate-13 lesson:

> **The structural coordinate has a scale.**

## 42 um support: still no natural length

For the tight `0.9` support, all tested objectives still run to a boundary.

```text
purity only        -> 5.0x boundary
0.1% noise         -> 5.0x boundary
1% noise           -> 0.2x boundary
3% noise           -> 0.2x boundary
```

At the 1% reference condition, even the full coherent edit has only

```text
phase shift              0.00226 rad
best-scalar residual      0.00650
```

So merely replacing `2 um` by `42 um` is not enough in this model.

## 356.6 um support: finite fixed points appear

At the `0.5` support, the structural coordinate spans about 29% of the input path and 54% of the integrated receipt field.

Now the grounded objective develops interior optima:

```text
observer floor     fixed scale     interior?
0.1%               1.5360x         yes
1.0%               0.6170x         yes
3.0%               0.2x            no
```

For the 1% reference observer:

```text
baseline support cable      356.60 um
fixed-point support cable   220.04 um
target gain                 1.297x baseline
target fraction             99.1365%
phase shift                  0.00985 rad
best-scalar residual         0.03048
```

The target becomes slightly *less pure* as a fraction, but its absolute transfer increases enough that the bounded observer's target/off-target-plus-noise score improves.

That is exactly why the observer grounding matters.

## 1098.8 um support: finite fixed points across all tested noise floors

The broad `0.1` support contains about 88% of the input route and 98% of its integrated receipt field.

All three tested grounded objectives now have interior fixed points:

```text
observer floor     fixed scale     target gain
0.1%               1.1746x         0.659x baseline
1.0%               0.8513x         1.425x baseline
3.0%               0.7248x         1.925x baseline
```

For the 1% reference case:

```text
baseline support cable      1098.79 um
fixed-point support cable    935.45 um
target gain                   1.425x baseline
target fraction              99.0784%
phase shift                    0.01198 rad
best-scalar residual           0.03545
```

Again, this is not “maximize resonance purity.” It is an observer-relative compromise between selectivity and a signal that still exists above the readout floor.

## The SighImageSuper control survives intact

Purity alone **never** finds an interior optimum for any of the three supports.

The broadest support is the strongest warning. At `5x` length:

```text
target fraction            99.999842%
target gain                 0.0000871 x baseline
```

So the operator can become almost perfectly pure by becoming almost perfectly silent.

That is the same conceptual trap exposed by `SighImageSuper`: the final surviving eigenmode can look maximally clean while absolute recoverable signal is disappearing.

Therefore:

```math
\text{mode dominance} \neq \text{useful information transmission}.
```

The bounded observer is not an afterthought. It changes what geometry counts as good.

## What Gate 14 actually adds

Gate 13 said:

> one sampled morphology interval is too microscopic a structural coordinate here.

Gate 14 adds:

> **finite grounded structural fixed points appear when the length variable acts coherently over enough of the transfer-relevant cable.**

In this one real morphology and one phenomenological membrane, the transition occurs somewhere between the tested `42 um` support and the tested `356.6 um` support for the 0.1% and 1% observer floors. The exact transition scale has not yet been identified and should not be turned into a biological constant.

The stronger candidate picture is now two-level:

```text
fast state iteration
    x_(t+1) = F_theta(x_t)
    -> preferred dynamical mode

slow structural iteration
    theta_(n+1) = G(theta_n, bounded consequence)
    -> supported geometry with no better neighbouring edit
```

A natural length is therefore not necessarily the length of one dendritic compartment. It can be a property of a **coherently changing structural support**.

## A correction before calling this credit assignment

The return field peaks at the driven terminal tip. That makes the next attack obvious.

We have shown that a transfer-derived spatial support can define a useful structural coordinate. We have **not** shown that multiplying forward occupancy by reciprocal return identifies the causal locus better than forward occupancy alone, distance from the input, or a generic distal growth zone.

The next control should compare equal-budget supports selected by:

```text
forward only
reverse only
forward x reverse
random contiguous region
matched-distance contiguous region
```

If the joint receipt wins after those controls, then “receipt chooses WHERE” earns more than metaphorical status. If not, the natural-length result can still stand while the receipt-credit story gets rejected.

## Claim boundary

Gate 14 establishes in this numerical system that:

- the real morphology supports a non-zero quasi-active transfer peak under the fixed phenomenological membrane;
- one approximately 42-um coherent support still has no tested interior grounded length optimum;
- a 356.6-um support has interior optima for two of three tested observer floors;
- a 1098.8-um support has interior optima for all three tested floors;
- purity alone always runs toward pathological isolation;
- coherent geometric edits produce materially non-scalar frequency-dependent operator changes.

It does not establish a biological growth unit, a universal natural dendritic length, a literal acoustic standing-wave condition, or a correct biological credit-assignment signal.
