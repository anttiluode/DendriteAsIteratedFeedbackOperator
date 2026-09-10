# Gate 9 — remove the privileged threshold

Gate 8 showed that Gate 7 was not a single numerical point, but the successful threshold interval was still narrow. Gate 9 therefore removes the threshold entirely.

The slow local write is now simply

```math
\Delta g = \eta \int V_{write}(t)^2\,dt.
```

Every event can write. At `eta = 30`:

```text
WA       0.03253
WB       0.05222
WAB      0.17526
Wsep     0.13864
WBsep    0.04947
```

So individual A and B effects are no longer artificially suppressed.

## Causal subtraction of the structural write

Use exact parallel-world inclusion/exclusion only as an analysis:

```math
\Delta g_{collision}
=
\Delta g_{AB}-\Delta g_A-\Delta g_B+\Delta g_0.
```

For simultaneous A+B:

```text
collision-specific write = 0.09051
```

For the temporally separated control, with its own time-matched B-alone world:

```text
collision-specific write = 0.05664
```

The matched nonlinear structural component is therefore `1.598x` the separated one.

That is much less spectacular than Gate 7's thresholded `12.26x`, but it is the more important result: the collision-specific persistent component survives after individual writes are explicitly allowed and subtracted.

## The later-response subtraction

Each world keeps its **own total physical write**. We do not create a magical world containing only the isolated counterfactual term.

Then all fast voltage and synaptic state is erased exactly. Every world receives the same passive unit-charge A probe. Only afterwards do we form the later-response interaction receipt:

```math
\chi_{AB}=R_{AB}-R_A-R_B+R_0.
```

At the B branch site:

```text
||chi matched||      0.019409
||chi separated||    0.012345
ratio                1.572x

peak |chi matched|   0.006236
peak |chi separated| 0.003974
ratio                1.569x
```

So a persistent later-routing component remains that is not explained by simply adding the two individual structural writes.

A write-scale attack from `eta = 10` through `100` keeps the B-site causal-response selectivity between `1.572x` and `1.618x`.

## Claim boundary

Inclusion/exclusion is a **counterfactual analysis across parallel worlds**. It is not a subtraction operation available to a biological neuron.

Gate 9 establishes a narrower statement:

> With a threshold-free local write in this toy, simultaneous nonlinear activity leaves a persistent structural and later-response interaction component that exceeds a time-separated control even after individual-event writes are explicitly accounted for.

The effect is modest, around `1.6x`, and should be treated as such.

The deterministic receipt is [`results/gate9_causal_write_receipt.json`](results/gate9_causal_write_receipt.json). Run:

```bash
python gate9_causal_write.py
```

The next serious move is the audited real morphology. If the same local-constraint -> global operator edit -> bounded visibility story survives there, the project is no longer resting on the geometry of a hand-built Y cable.
