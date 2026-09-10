# Gate 13 — the natural-length rule hits a discretization-scale wall

Gate 12 found finite, observer-grounded length fixed points in the small Y-cable toy. Gate 13 transfers the same question to the pinned human L2/3 pyramidal morphology already used by `Operaattori`.

The result is **not** a clean replication.

That is useful.

## What is real, and what is still a toy?

The morphology is the pinned cell-1125 ASC from `ido4848/FCI` at commit `75ad8b4d81a7f51bf888b30650c543592340db06`. After retaining soma and dendritic sections, the compiled point tree has:

```text
12,514 dendritic nodes including the soma root
12,513 dendritic edges
20,742.41 um total dendritic cable
99 branch nodes
```

The membrane remains a uniform phenomenological quasi-active model. This is therefore a **real-geometry operator test**, not a fitted human-neuron simulation.

The longest soma-to-tip cable path in this tree is about

```text
1,242.17 um across 750 sampled path nodes.
```

The baseline operator has an interior non-zero transfer peak at

```text
omega* = 175.869 rad/s
```

and, for the same five-tone challenge used in Gate 12, the target component already occupies

```text
99.3765% of AIS output power.
```

So the real-geometry test begins from an operator that is already extremely selective under the chosen uniform quasi-active membrane.

## WHERE comes from forward × reciprocal return

The candidate edit is **not** chosen by looking for the best length derivative.

At `omega*`, Gate 13 computes:

```text
forward = response from the selected distal tip
return  = reciprocal small-signal response from the AIS
where_j = |forward_j| |return_j|
```

restricted to the actual soma-to-input path.

The maximum internal overlap picks node `10531`, parent `10530`.

Its physical segment is only

```text
1.98457 um long
radius = 0.745 um
```

or about **0.160% of the 1,242 um input path**.

That number turned out to matter.

## The Gate-12 fixed point does not survive at one sampled segment

We swept that one real segment over

```text
0.2x ... 5x its measured length
```

while keeping radius and membrane densities fixed. Length therefore changes local membrane area, capacitance, leak and quasi-active conductance together, while axial conductance changes inversely with length.

### Purity-only objective

Purity walks all the way to the upper boundary:

```text
scale                         5.0x
target fraction               99.3867%
target gain vs baseline       0.9853x
interior fixed point?         no
```

### Grounded bounded-observer objective

We then used the Gate-12 score

```math
Q_\sigma(L)=\frac{P_{target}(L)}{P_{off}(L)+\sigma^2}.
```

The result is still boundary-seeking rather than fixed-point seeking.

```text
noise floor        preferred boundary
0.1%               5.0x
0.3%               5.0x
1.0%               0.2x
3.0%               0.2x
10%                0.2x
```

There is **no finite interior natural length** in this tested one-segment coordinate.

That is the Gate-13 result.

## Why this is not evidence against geometry-as-operator

The more informative measurement is how little this edit changes the global transfer.

At the 1% reference condition, shrinking the chosen segment from `1.985 um` to `0.397 um` produces only

```text
target-transfer phase shift          1.075e-4 rad
five-tone best-scalar residual       3.047e-4
```

In Gate 12, the corresponding residual was about `2.55e-2`.

So the real-morphology transfer exposed a scale mismatch in our question:

> **a MorphIO sampling interval is not automatically the biological analogue of Berglund's channel length.**

Gate 12 let one parameter represent a substantial fraction of a branch. Gate 13 let one parameter represent roughly two micrometres out of a 1.24 mm path. Even a five-fold edit changes the total path length by less than one percent.

The failure therefore says something more precise than “natural length is wrong”:

> **The natural-length hypothesis is scale-sensitive. A coordinate defined by file discretization can be too microscopic to expose a meaningful operator fixed point.**

## The next attack

Do not tune the membrane until a fixed point appears.

Change the structural coordinate.

The Berglund analogue is a coherent channel, and biological structural plasticity can alter spines, branch stretches, terminal growth zones and whole arbor segments. The next gate should therefore promote the receipt from a point to a **spatial support**:

```text
receipt field
    -> choose a contiguous supported cable region
    -> coherently lengthen / shorten that region
    -> ask the same bounded observer again
```

This preserves the useful Gate-11 decomposition:

```text
receipt chooses WHERE
small structural perturbation discovers WHICH WAY
consequence decides WHETHER TO KEEP IT
```

but stops pretending that one numerical sample interval is necessarily the natural unit of growth.

## Claim boundary

Gate 13 establishes only that:

- the pinned real dendritic morphology can be compiled into the same frequency-domain operator experiment;
- the chosen forward/return overlap identifies a real physical path segment without consulting the length objective;
- changing that one approximately 2-um segment produces only a tiny non-scalar operator deformation;
- neither purity nor the tested bounded-observer scores produce an interior fixed point over `0.2x ... 5x` for that coordinate.

It does **not** establish that real dendrites seek eigenmodes, that bAPs provide this receipt, that structural plasticity follows this optimizer, or that dendrites grow to literal acoustic wavelengths.

The next falsifiable question is whether a receipt-sized **coherent cable region**, rather than one morphology-file interval, has a finite grounded structural fixed point.
