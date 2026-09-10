# Gate 15 — Equal-budget support response

## Question

Gate 14 looked as if an extended dendritic support might define a functional coherence scale: a ~42 um receipt support had no interior natural-length fixed point, while ~357 um and ~1099 um supports did.

But Gate 14 changed every segment in a support by the same *fractional* amount. A larger support therefore received a larger total structural edit. Gate 15 asks the harder question:

> **If every support gets exactly the same total cable-length change budget, does operator leverage still switch on only above a characteristic spatial width?**

Define

```math
\Gamma_B(W)=
\max_{\Delta L:\,\sum_j |\Delta L_j|=B}
\left[Q_\sigma(\theta+\Delta L)-Q_\sigma(\theta)\right],
```

with the unedited state included so a support is allowed to say `stay`.

For this gate:

```text
fixed structural budget B = 4 um total |Delta L|
observer floor             = 1% of baseline target power
support widths             = 20, 40, 80, 160, 320, 640, 1000 um
real morphology            = pinned human L2/3 cell 1125
baseline carrier           = 175.868675 rad/s
baseline target fraction   = 99.3765%
```

The full deterministic output is uploaded by GitHub Actions. A compact frozen receipt is [`results/gate15_support_response_summary.json`](results/gate15_support_response_summary.json).

## Result 1 — the proposed sharp support threshold does **not** survive the equal-budget control

For the joint `|forward| x |return|` selector, relative bounded-score improvement is:

```text
requested W     actual cable      Gamma_B / Q0
20 um           21.35 um          +0.5762%
40 um           41.29 um          +0.5759%
80 um           81.18 um          +0.5758%
160 um         161.94 um          +0.5802%
320 um         320.92 um          +0.5954%
640 um         640.06 um          +0.6508%
1000 um        999.89 um          +0.6828%
```

There is no turn-on between the ~42 um and ~357 um Gate-14 supports. Useful local leverage is already present at the smallest width tested here, and the curve changes smoothly.

The descriptive half-max onset is therefore the **smallest tested support** at every attacked frequency.

This means Gate 14 established something narrower than the first interpretation:

> **Extended coherent geometry can strongly alter the operator under the Gate-14 scaling protocol, but Gate 14 did not establish a minimum operator-defined spatial support scale.**

The apparent Gate-14 scale transition was at least partly confounded by edit magnitude growing with support size.

## Result 2 — `forward x return` has not earned credit-assignment status

Five support selectors received the same physical width and the same 4 um structural budget:

- forward-only;
- reverse-only;
- joint `|forward| x |return|`;
- deterministic random contiguous support;
- a gross-distance-matched, low-joint control.

At the baseline frequency, joint and forward select the same support exactly at most widths and nearly the same at the rest. More damaging to the strong receipt story, controls often do better.

Examples:

```text
W = 20 um
joint       +0.5762%
random      +1.0648%

W = 160 um
joint       +0.5802%
reverse     +0.7687%
random      +0.7747%

W = 320 um
joint       +0.5954%
reverse     +0.7979%
random      +0.8699%

W = 1000 um
joint       +0.6828%
reverse     +0.7335%
distance    +0.7337%
```

So the return field is not currently identifying a uniquely useful growth support. On this path, the joint field is largely dominated by the forward field near the driven route.

That is a useful negative result: the Gate-11 receipt can mark an event-conditioned return, but **the present `forward x return` support heuristic is not yet a demonstrated structural credit signal.**

## Result 3 — every winning local edit is shortening

Across every width and every selector at the baseline carrier, the locally best equal-budget move is `shorter`.

The same is true in the joint support sweep at `0.5x`, `1x`, and `2x` the baseline carrier.

That tells us not to over-read the local response as a resonance-matching growth rule. Around the present morphology, the bounded objective has a broadly shared local direction: reduce path metric slightly.

Gate 14's interior optima can still exist farther along that direction; Gate 15 says the *local derivative-like signal* is not specific to a large coherent support.

## Result 4 — frequency changes the size of the effect, not the onset

Joint-selector equal-budget curves:

```text
0.5 x omega*: 20 um +0.6406%, best ~640 um +0.7178%
1.0 x omega*: 20 um +0.5762%, best ~1000 um +0.6828%
2.0 x omega*: 20 um +0.2327%, best ~640 um +0.2984%
```

A local cylindrical quasi-active cable coordinate was also computed from

```math
\gamma^2(\omega)=r_a\,y_m(\omega),
```

and every support is reported both in micrometres and accumulated `|gamma dx|`.

But because the half-max onset is already at the first physical width for all three frequencies, there is no threshold to collapse into a common electrotonic coordinate. The electrotonic-scale hypothesis remains open, not supported by this gate.

## What Gate 15 kills

Gate 15 does **not** support the strong statement:

> `the operator defines a minimum spatial scale below which geometry edits do nothing.`

Not over the tested 20–1000 um range with this structural budget and this model.

It also does not support:

> `forward x return uniquely selects the structurally useful region.`

## What survives

Several earlier results remain intact:

1. changing dendritic length changes the complex frequency-dependent operator non-scalarly;
2. bounded-observer grounding prevents the purity-through-silence failure;
3. real morphology admits observer-grounded interior geometry optima for some extended coordinates;
4. the useful structural coordinate need not be a morphology-file sample interval;
5. the support-response function `Gamma_B(W)` is a clean measurable object, even when it falsifies the threshold story.

The important correction is therefore:

```text
Gate 13: one sampled point was a bad structural coordinate.
Gate 14: larger coherent supports can produce strong finite optima.
Gate 15: after equalizing structural budget, there is no sharp support turn-on.
```

## Next attack

The smallest Gate-15 width is still ~20 um, whereas Gate 13's microscopic coordinate was ~2 um. The next clean experiment is not to rescue a threshold by tuning parameters. It is to bridge that decade directly with a much smaller fixed structural budget:

```text
2, 5, 10, 20, 40, ... um
```

and repeat after morphology re-sampling.

If leverage simply continues smoothly to the single-edge limit, the coherence-threshold idea should be retired.

If a reproducible onset survives **equal budget + resampling**, then and only then does a functional spatial scale become an earned object.

## Claim boundary

The morphology is real. The membrane kinetics, target tones, bounded observer and structural search are phenomenological. `Gamma_B(W)` is an operator-sensitivity experiment, not a demonstrated biological growth law. The local electrotonic coordinate is a diagnostic from cylindrical cable theory, not a literal dendritic wavelength.

The point of the gate is exactly that the attractive interpretation was allowed to lose.
