# Dendrite As Iterated Feedback Operator

This repo asks a stronger question than whether a passive dendrite is a leaky cable:

> **Can dendritic structure + membrane kinetics act as an eigenmode purifier — rejecting most of a mixed temporal input while making one branch-selected non-zero mode dominate the bounded output?**

In the current numerical toy, **yes in the narrow filtering sense.**

The repo then attacks the stronger follow-up: does that purified mode specifically control a local nonlinear event? The first nonlinear gate works, but an ablation shows that its selectivity is **not yet causally attributable to resonance**. That failure is preserved.

The working ontology is:

```text
local geometry + membrane kinetics
        ->
frequency-dependent propagation operator
        ->
mode selection / purification
        ->
local nonlinear event
        ->
state-dependent operator
        ->
slow material edit
        ->
future propagation changes
```

The first four lines are now executable in pieces. The causal link between **purification** and **nonlinearity** is the next target.

---

# Gate 5 — the eigenmode purifier

The passive cable is the control. It obeys

```math
C\dot V=-GV+I
```

and its autonomous modes simply decay. In the passive toy, the distal-to-soma transfer peaks at DC. Repeated evolution eventually leaves the slowest decaying direction: useful spectral selection, but basically diffusive forgetting.

Gate 5 adds one phenomenological restorative state per dendritic compartment:

```math
C\dot V=-GV-g_qz+I,
```

```math
\tau_q\dot z=V-z.
```

Under sinusoidal drive the local quasi-active contribution is

```math
Y_q(\omega)=\frac{g_q}{1+i\omega\tau_q},
```

so the global transfer becomes

```math
H(\omega)=
\left[
G+i\omega C+
\operatorname{diag}Y_q(\omega)
\right]^{-1}.
```

This is a phenomenological restorative membrane model, **not** a fitted HCN model and not a claim that a biological dendrite is an acoustic waveguide.

## Two branches, two non-zero pass bands

GitHub Actions gives:

```text
branch A
    peak omega       0.125557
    peak / DC gain   29.91x
    Q-like           1.654

branch B
    peak omega       0.068077
    peak / DC gain   34.31x
    Q-like           1.388
```

The full `V + z` system is stable and contains `12` complex poles. The least-damped measured pair is approximately

```text
-0.03997 +/- i 0.05469
```

So this is categorically different from the passive RC control's purely real decays, while remaining a dissipative electrical state-space model.

## The actual purification test

Feed five equal-energy temporal components into a distal tip. The target resonant component starts with only

```text
20% of input energy.
```

At the soma:

```text
branch A active target fraction   97.61%
branch B active target fraction   96.03%

branch A passive target fraction   0.64%
branch B passive target fraction   4.36%
```

So the toy really does this:

```text
mixed temporal input
        ->
distributed dendritic operator
        ->
most components rejected
        ->
selected non-zero mode dominates output
```

That is the precise operational meaning of **eigenmode purifier** used here.

The frozen receipt is [`results/active_receipt.json`](results/active_receipt.json).

---

# Gate 6 — can the selected activity enter a nonlinear regime?

A bounded generic regenerative hotspot was placed on branch A. Two same-branch sources were phase/amplitude calibrated at branch A's preferred frequency so that the matched pair reaches the hotspot coherently.

At the default operating point:

```text
matched linear hotspot peak      0.0800
quadrature hotspot peak          0.05657
single hotspot peak              0.0400
nonlinear knee                   0.0600
```

The soma nonlinear residual is:

| condition | residual RMS |
|---|---:|
| matched | `4.2039e-3` |
| quadrature | `2.8811e-3` |
| high off-band | `6.93e-10` |
| spatially separated | `1.60e-6` |
| single | `1.48e-6` |

So coherent same-branch activity can indeed be converted into a selective nonlinear event. The hardest control is quadrature, and matched beats it by `1.459x`.

This hotspot is deliberately generic. **It is not an NMDA kinetic model.**

The frozen receipt is [`results/nonlinear_receipt.json`](results/nonlinear_receipt.json).

---

# Gate 6b — attack the result rather than celebrate it

Two obvious attacks were run.

## Attack 1: was the result only the chosen threshold?

Move the nonlinear knee while freezing carrier, phases, amplitudes, substrate and geometry:

```text
knee 0.045 -> matched/quadrature   1.151x
knee 0.050 ->                      1.197x
knee 0.055 ->                      1.278x
knee 0.060 ->                      1.459x
knee 0.065 ->                     41.822x
knee 0.070 ->                    261.534x
knee 0.075 ->                   1317.572x
```

The matched advantage therefore is not confined to one hand-picked threshold.

## Attack 2: is resonance actually necessary?

Remove all quasi-active conductance while keeping the **identical raw Gate-6 protocol**: same carrier, source amplitudes, phase, addresses, hotspot and nonlinear knee. No passive re-tuning.

The passive ablation still gives:

```text
matched hotspot peak             0.06633
matched nonlinear residual       0.002927
quadrature residual              1.93e-8
```

The active resonant substrate makes the matched hotspot `1.206x` larger and the matched nonlinear residual `1.436x` larger.

But the passive substrate still performs strong coherent coincidence selection.

Therefore:

> **Gate 6 does not yet establish that resonance caused the nonlinear selectivity.**

That matters. A clean single carrier can add coherently even in a passive cable. So the next experiment must make the purifier indispensable rather than merely present.

Full attacker write-up: [`GATE6.md`](GATE6.md). Frozen receipt: [`results/gate6b_attacker.json`](results/gate6b_attacker.json).

---

# The older structural results still matter

## Passive Sigh-style iteration

A distal pulse begins with modal effective dimension `10.005` and falls under silent passive evolution to `1.019`; the late normalized state reaches cosine `0.99863` with the slowest mode.

That is the SighImageSuper bridge, but only as a passive control.

## Branch != scalar weight

Equal charge at two distal tips produces soma kernels peaking at steps `79` and `93`. Even after fitting the best scalar map, `18.23%` of one waveform remains unexplained.

A branch carries a temporal transfer function, not merely a multiplier.

## One local edit -> global low-rank operator edit

For one axial conductance edit

```math
G' = G + \delta g\,bb^T
```

and

```math
H(s)=[G+sC]^{-1},
```

Sherman-Morrison gives

```math
\Delta H
=-\frac{\delta g Hbb^TH}
        {1+\delta g b^THb}.
```

The passive numerical update has effective rank exactly `1.0` to roundoff.

At fixed frequency, one local quasi-active conductance edit also gives a global transfer edit with effective rank

```text
1.0000000000009
```

and `s2/s1 ~ 1e-14`.

So the recurring structural claim survives both passive and quasi-active versions:

> **A local physical change can make a distributed but low-dimensional edit to the global response operator.**

---

# Why this connects Sigh, the resonant cavities, and the neuron work

`SighImageSuper` showed:

> **the operator determines which distinctions persist.**

The resonant-cavity work showed:

> **sparse local constraints can compile a dense global Green's function.**

This repo now adds:

> **a dendrite-like distributed electrical operator can strongly enrich a non-zero temporal mode before a bounded readout.**

The important claim is not "dendrites are little acoustic cavities."

It is:

```text
weights / parameters = local physical constraints
operator             = propagation implied by those constraints
state                = activity occupying its modes
computation          = which modes reach which nonlinear regions/readouts
learning             = local changes that alter future propagation
```

That is much closer to the original intuition than flattening a dendrite into one scalar `w_ij`.

---

# Next gate — make purification causally necessary

Do **not** give the nonlinear hotspot a clean carrier that the passive cable can also sum.

Reuse Gate 5's five equal-energy temporal components, so the target mode starts at only `20%` of the input. Give the active and passive substrates the exact same broadband waveform and the exact same hotspot.

Then ask whether:

```text
ACTIVE
broadband mixture
    -> branch purifies target mode
    -> target crosses nonlinear regime
    -> bounded nonlinear response

PASSIVE ABLATION
same broadband mixture
    -> no target purification
    -> no comparable selective nonlinear event
```

No re-tuning after ablation.

If that survives, then we can finally make the stronger sentence:

> **The dendritic operator purified a mode, and that purified mode triggered the nonlinear event.**

After that:

```text
selected mode
    -> nonlinear state change
    -> Jacobian changes
    -> slow local write
    -> future resolvent changes
```

and then transfer the experiment to the audited real `Operaattori` morphology.

---

# Run

```bash
python -m pip install -e .[dev]
pytest -q
python experiment.py
python active_experiment.py
python nonlinear_experiment.py
python gate6b_attacker.py
```

CI runs the full sequence on Python 3.10 and 3.12.

See [`THEORY.md`](THEORY.md), [`GATE6.md`](GATE6.md), and the frozen `results/` receipts.

---

## Claim boundary

Established in this numerical toy:

- passive cable evolution performs diffusive mode selection;
- branch transfer is not reducible to one scalar weight;
- one local edge/conductance edit can make a global rank-one resolvent change;
- a restorative quasi-active extension gives stable complex poles and branch-specific non-zero pass bands;
- a deliberately tuned quasi-active branch can enrich one of five equal-energy temporal components from `20%` input energy to about `96-98%` of bounded-output energy;
- coherent same-branch activity can drive a local regenerative nonlinearity;
- the Gate-6 matched advantage survives a nonlinear-threshold sweep.

Not established:

- that resonance is necessary for the current two-tone nonlinear coincidence result — the passive ablation disproves that claim for Gate 6;
- that real dendrites routinely achieve the toy's purification strength;
- that HCN alone matches the toy parameters;
- that biological dendrites are acoustic/RF waveguides;
- that NMDA automatically implements resonant mode competition;
- general learning, transformer equivalence, energy advantage, or brain-level computation.

Those are experiments, not conclusions.
