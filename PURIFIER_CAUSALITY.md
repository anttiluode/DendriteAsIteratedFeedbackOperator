# Purifier causality audit

The central thought in this repository is stronger than "a passive dendrite forgets high spatial modes faster than low ones."

The claim worth testing is:

> **A distributed dendritic operator can progressively enrich a particular non-zero temporal mode while rejecting other components, so downstream nonlinear machinery sees a cleaner dynamical object than the distal input contained.**

Gate 5 already showed this at the soma. This audit asks whether the purification is actually happening **along the dendrite**, and then separates that fact from the later nonlinear/plasticity gates.

## 1. Passive decay is only the control

For a passive cable,

```math
C\dot V=-GV,
```

all autonomous modes decay. Repeated propagation eventually leaves the slowest decaying direction. That is useful mathematics, but it is fundamentally a diffusive/low-pass story.

The stronger object begins only after membrane dynamics make the transfer genuinely frequency selective. In the quasi-active Gate-5 toy,

```math
C\dot V=-GV-g_qz+I,
```

```math
\tau_q\dot z=V-z,
```

which gives

```math
H(\omega)=
\left[G+i\omega C+\operatorname{diag}\frac{g_q}{1+i\omega\tau_q}\right]^{-1}.
```

The resulting state operator has stable complex poles and a non-zero pass band. That is the first point at which **eigenmode purifier** is an appropriate operational phrase in this repo.

## 2. The purifier is distributed, not a soma artifact

A distal input is sent into branch A. At each successive branch compartment we measure five non-DC temporal components:

```text
0.35 w_peak, 0.60 w_peak, w_peak, 2.5 w_peak, 6.0 w_peak
```

They have equal input amplitude, so the target starts as one component among five. No nonlinear threshold or hotspot is used to select the measurement location.

The striking result is that the target fraction rises as the signal propagates proximally through the **active** branch:

| branch site | active target fraction | passive target fraction | active/passive purity |
|---:|---:|---:|---:|
| distal tip 12 | 43.97% | 16.81% | 2.62x |
| 10 | 55.93% | 13.68% | 4.09x |
| 8 | 66.27% | 10.73% | 6.18x |
| 6 | 74.23% | 8.22% | 9.03x |
| 4 | 80.46% | 6.36% | 12.66x |
| proximal 1 | **83.74%** | **3.91%** | **21.41x** |

So the phenomenon is not:

```text
mixed signal -> arbitrary soma readout happens to like one frequency
```

It is much closer to:

```text
mixed distal signal
    -> repeated propagation through distributed membrane dynamics
    -> off-target components progressively suppressed relative to target
    -> cleaner target-dominated state near the branch exit
```

That is the SighImageSuper intuition in a more literal dendritic form: **the route itself changes the modal composition of the state.**

The deterministic receipt is [`results/purifier_hotspot_receipt.json`](results/purifier_hotspot_receipt.json).

## 3. Purification is not the same thing as amplification

There is an important counter-result.

At every sampled site, the matched passive cable has **more total broadband gain power** than the quasi-active branch. For example, at proximal site 1:

```text
active total gain power   0.00194
passive total gain power  0.02677
```

Yet the active branch has `83.74%` of that output power concentrated in the target component, while the passive branch has only `3.91%`.

So the right picture is not an active branch indiscriminately boosting its preferred signal.

It is:

> **the active branch pays amplitude to buy selectivity.**

That matters for the next nonlinear experiment. A memoryless voltage threshold could easily fire more strongly in the passive cable simply because the passive waveform is larger overall, even though it is spectrally dirtier.

Therefore a claim like

```text
resonance -> larger voltage -> NMDA -> learning
```

would not isolate the purifier mechanism.

## 4. The nonlinear/plasticity gates are real but currently parallel

The concurrent Gates 6–9 on `main` establish a different sequence in the same Y-cable family:

```text
local voltage-dependent nonlinearity
    -> state-dependent Jacobian
    -> local persistent conductance write
    -> exact fast-state erasure
    -> changed later route
```

Those are useful results. Gate 6 shows that the fast nonlinear state changes the instantaneous transfer operator; Gate 9 leaves a modest but persistent collision-specific later-response component after individual writes are explicitly accounted for.

But those experiments do **not yet prove** that Gate-5 modal purification is what caused the nonlinear event or the later write.

A separate two-tone attacker made this especially clear: when the quasi-active resonance was removed but the same coherent two-source protocol was retained, a passive cable could still sum the clean carrier and drive a local threshold strongly. Coherent coincidence and resonant purification can coexist without the latter being causally necessary.

The project should therefore keep two statements separate for now:

```text
A. dendritic distributed mode purification        established in toy
B. nonlinear state -> operator change -> write   established in toy
A causes B                                        not established yet
```

## 5. The next decisive experiment

The next gate should make modal composition, rather than raw amplitude, the controlled variable.

Use a genuinely multimode distal input. Measure at a predeclared proximal branch location where the linear diagnostic already shows strong active/passive purity contrast. Then compare active and passive substrates under **two budgets**:

1. **equal distal input:** the natural physical comparison;
2. **equal local RMS/power:** rescale only for an attacker so the nonlinear element cannot win merely because one substrate delivers more total voltage.

At equal local power the key contrast becomes:

```text
same total local activity
        but
active: target-dominated temporal composition
passive: broadband / low-frequency-dominated composition
```

Now give both worlds the **same dynamical nonlinear detector**. It must have temporal state; a memoryless `f(V)` cannot distinguish two signals solely because their spectra differ once instantaneous amplitude statistics are matched.

The cleanest toy detector would be a local conductance/gating system with its own finite time constant and voltage dependence, followed by the same structural-write test already used in Gates 7–9. The detector parameters must be frozen before the active/passive comparison and must not be separately tuned to each substrate.

If the active branch then produces a larger nonlinear state or persistent write at equal local power, while the passive branch does not, the causal sentence becomes much stronger:

> **The dendritic operator purified a mode, and the purified temporal structure changed what the nonlinear branch machinery did.**

## 6. What the idea is becoming

The most interesting version is no longer "dendrites are resonant waveguides."

It is more general and more testable:

```text
morphology + membrane kinetics
        -> distributed operator
        -> progressive modal purification
        -> local state-dependent nonlinearity
        -> local constraint edit
        -> global transfer-function edit
        -> different future purification / routing
```

That is a closed computational loop.

And it puts the old Sigh sentence into dendritic language:

> **A dendrite does not merely attenuate a signal on its way to the soma. Its distributed dynamics can decide which temporal distinctions are allowed to survive the trip.**

The next serious endpoint remains the audited real `Operaattori` morphology. The Y cable is useful precisely because it lets each causal link be attacked before the same questions are asked of a reconstructed cell.
