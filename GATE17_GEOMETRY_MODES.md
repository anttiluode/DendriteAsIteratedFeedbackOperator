# Gate 17 — Geometry eigenspectrum and receipt projection

Gates 15-16 killed the proposed sharp minimum support width. Geometry affects the transfer operator all the way down to the smallest tested physical support, so the remaining question is not **where sensitivity turns on**. It is whether many locally sensitive geometric coordinates organize into a smaller set of consequential collective directions.

Gate 17 tests that directly on the same pinned human L2/3 morphology used in Gates 13-16.

## Protocol

The selected soma-to-distal-tip route is `1242.169 um` long in the `12,514`-node morphology. Segment lengths on that path are represented in a length-weighted orthonormal cosine basis in **log length**:

```math
\delta \log L(x)=\sum_k a_k v_k(x),
```

with

```math
\langle v_i,v_j\rangle_L
=\int v_i(x)v_j(x)\,\frac{dx}{L_{\rm path}}
=\delta_{ij}.
```

Thus the same coefficient amplitude gives the same length-weighted RMS fractional geometry perturbation regardless of how many morphology samples happen to represent a region.

The gate measures central finite-difference derivatives for `16` geometry modes at amplitude `eps = 0.01`, then builds the full `8 x 8` Hessian in the first eight smooth modes:

```math
g_k=\frac{\partial Q_\sigma}{\partial a_k},
\qquad
K_{ij}=\frac{\partial^2 Q_\sigma}{\partial a_i\partial a_j}.
```

Because the objective is **maximized**, negative Hessian eigenvalues are locally concave directions. Positive curvature would indicate a local minimum along that direction, not a stable maximizing direction.

The same gate also attacks the bra-ket/receipt proposal. A fixed phenomenological post-event conductance state — a proximal `120 um` conductance change plus an AIS shunt — produces

```math
|\delta r\rangle
=|r_{\rm post}\rangle-|r_{\rm pre}\rangle.
```

The following real spatial densities are projected onto exactly the same geometry modes and compared against the **measured** geometry gradient:

```text
Re <delta r|e> locally     event phase overlap
|delta r| |e|              event magnitude product
|e|^2                       forward-only control
|delta r|^2                 return-only control
Re <r_pre|e> locally       reciprocal-return control
```

Circular shifts of each density along the path provide a deterministic spatial-null test.

The full GitHub Actions outputs are archived by CI. The compact frozen receipt is [`results/gate17_geometry_modes_summary.json`](results/gate17_geometry_modes_summary.json).

## Result 1 — the geometry gradient is strongly low-order, but mostly global shortening

The first four smooth cosine modes contain

```text
97.8809% of total gradient energy.
```

The last four fine modes contain only

```text
0.2627%.
```

That initially looks like strong low-dimensional structural organization. But the attack matters: **mode 0 alone contains 96.7318% of all gradient energy.** Its relative derivative is

```text
mode 0: -2.21070 per unit mode amplitude
```

whereas the remaining modes are much smaller. This is consistent with the broad `shorter` direction already found in Gates 15-16.

Among the *non-uniform* modes, the first eight still contain about `70.40%` of the remaining gradient energy, so there is residual large-scale structure. But Gate 17 does not justify claiming that the 16-mode gradient is a rich low-dimensional shape code: most first-order leverage is a common path-scale direction.

## Result 2 — the smooth geometry subspace is locally concave

All eight eigenvalues of the measured `8 x 8` Hessian are negative:

```text
relative curvature eigenvalues, ordered by |lambda|

-15.60789
 -2.69907
 -2.39069
 -2.31634
 -2.22050
 -1.97849
 -1.89355
 -1.85725
```

The leading curvature direction is about `5.78x` larger in magnitude than the second and accounts for about `87.67%` of squared curvature energy. It is almost the uniform geometry mode (`0.9863` coefficient on mode 0).

The complete eight-mode curvature is **not rank one**:

```text
effective rank of |lambda| spectrum   5.23
95% Frobenius-energy rank             5 / 8
negative eigenvalues                  8 / 8
positive eigenvalues                  0 / 8
```

So the useful statement is narrower:

> **The tested smooth geometry subspace is locally concave, with one dominant almost-global curvature direction plus several weaker collective shape directions.**

This is already different from Gate 16's pointwise sensitivity result. Local sensitivity can exist everywhere while curvature organizes collectively in geometry-mode space.

## Result 3 — a collective eight-mode stationary candidate survives direct nonlinear replay

The Hessian and gradient predict the local Newton candidate

```math
|a_N\rangle=-K^{-1}|g\rangle.
```

Its coefficient vector is

```text
[-0.15544, -0.05370, +0.09072, +0.10568,
 +0.01336, +0.03552, +0.02222, +0.01213]
```

with `L2 = 0.22029` and maximum absolute coefficient `0.15544`.

This candidate was **not accepted from the quadratic approximation alone**. A counted nonlinear replay scanned `0, 0.25, ..., 1.5` times the candidate:

```text
fraction of Newton step       relative Q change
0.00                           0.0000%
0.25                          +8.0580%
0.50                         +13.5497%
0.75                         +16.5383%
1.00                         +17.3321%   <-- best counted point
1.25                         +16.3631%
1.50                         +14.0849%
```

The directly measured first-gradient norm in the eight-mode subspace then falls from

```text
2.23436  ->  0.43032
```

or to

```text
19.26% of its baseline value.
```

So the collective candidate is not merely a Hessian artifact. It lies near a real interior improvement peak along the predicted collective direction and is substantially more stationary than the original geometry within this eight-mode subspace.

This is the strongest surviving version of the **natural geometry** idea so far:

> **Not every segment has its own natural length. A collection of locally sensitive length coordinates can possess a finite observer-grounded optimum in a collective geometry subspace.**

The result is still dominated by overall shortening. Gate 18 must therefore hold global path scale fixed and ask whether a genuine *shape-only* collective fixed point survives.

## Result 4 — the bra-ket receipt is a useful decomposition, but this receipt fails credit assignment

The global normalized overlap between the event-conditioned return residual and forward eligibility is only

```text
|<delta r|e>|^2 / (<delta r|delta r><e|e>) = 0.02843.
```

The post-event return residual is also extremely close to a scalar deformation of the ordinary reciprocal return:

```text
non-scalar residual = 0.00775
```

so this simple post-event conductance state adds less than one percent of a new spatial return direction.

The event phase-overlap projection nevertheless has a seemingly attractive raw correlation with the measured geometry gradient:

```text
corr(mode projection, measured gradient) = +0.7743
```

but the spatial-null attack kills the strong interpretation:

```text
circular-shift null median |corr| = 0.8174
empirical p                         = 0.9531
```

A shifted version of the receipt is typically **at least as predictive**. The reciprocal-return control gives essentially the same magnitude (`|corr| = 0.7733`). The event-conditioned overlap therefore has not earned structural credit-assignment status.

This is exactly why the Dirac notation is helpful without being magical:

```math
\langle r|e\rangle
```

is a well-defined overlap, but overlap is not automatically causal permission; and

```math
|e\rangle\langle r|
```

is a well-defined rank-one operator in activity space, but it is not automatically the physical map from activity to dendritic geometry.

The still-missing arrow is

```math
(|e\rangle,|r\rangle)
\longrightarrow
|\delta \ell\rangle
\longrightarrow
\Delta H.
```

Gate 17 measures the first arrow against the actual geometry derivative and, for this receipt construction, it fails the spatial-specificity control.

## What Gate 17 establishes

Within one real morphology and one phenomenological quasi-active/observer model:

- smooth log-length geometry modes provide a clean, discretization-independent structural coordinate;
- first-order geometry sensitivity is heavily concentrated in a global shortening mode;
- the first eight smooth geometry modes have an all-negative local Hessian;
- one curvature mode dominates, but the full curvature has effective rank about `5.23` rather than one;
- the eight-mode Newton direction gives a **+17.33%** directly replayed bounded-score improvement and reduces the local gradient norm by about `80.7%`;
- the tested event-conditioned receipt overlap does **not** beat spatial-shift or reciprocal-return controls.

## What Gate 17 does not establish

It does not show that:

- dendrites are quantum systems because bra-ket notation is convenient;
- real neurons compute a Hessian, Newton step, or cosine geometry basis;
- the actual biological morphology was produced by optimizing this objective;
- the geometry-change operator is globally low rank;
- the event-conditioned return is a useful biological credit signal;
- the non-uniform shape degrees of freedom matter independently of global shortening.

## Next attack — remove the easy scale direction

Gate 18 should remove mode 0 completely and constrain every allowed perturbation to preserve global path scale to first order:

```math
\langle v,1\rangle_L=0.
```

Then repeat the gradient/Hessian/validation experiment in shape-only coordinates.

If the interior collective optimum disappears, Gate 17 was mainly another expression of the broad shortening direction already visible in Gates 15-16.

If a finite optimum survives in zero-mean geometry modes — especially across several paths or cells — then the stronger object becomes earned:

> **Geometry has collective shape modes whose equilibrium changes the operator even after global size is held fixed.**

That is the clean next bridge to the SighImageSuper idea of many microscopic coordinates collapsing onto a smaller set of consequential modes.

## Claim boundary

The morphology is real. The membrane kinetics, observer objective, post-event conductance state and chosen cosine geometry basis are phenomenological. The finite-difference Hessian is a local property of that model, not a biological learning algorithm. The receipt test is intentionally negative where its controls win.

The useful result is the separation:

```text
local geometry sensitivity
        !=
collective structural curvature
        !=
credit assignment.
```

Gate 17 finds evidence for the second while refusing to infer the third.
