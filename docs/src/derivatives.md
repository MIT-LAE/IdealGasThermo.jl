# [Thermodynamic derivatives](@id derivatives)

The pure core is generic over `Real`, so
[ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) differentiates through
every property and process out of the box. This page explains *how* the
derivatives are propagated — the closed-form property rules, what a "partial"
is, and the implicit-function-theorem rules that differentiate the Newton
inversions **without re-running the solver**.

## Closed-form property derivatives

The forward properties have exact derivatives:

```math
\frac{dh}{dT} = c_p, \qquad
\frac{ds^0}{dT} = \frac{c_p}{T}, \qquad
\frac{dc_p}{dT} = \Ru\!\left(-2a_1T^{-3} - a_2 T^{-2}
+ a_4 + 2a_5T + 3a_6T^2 + 4a_7T^3\right).
```

## Forward-mode AD and "partials"

Forward-mode AD evaluates a function on a **dual number** — a value carried
together with its derivatives. A dual with ``N`` *partials* is

```math
x + \sum_{i=1}^{N} \dot{x}_i\,\epsilon_i,
```

where each partial ``\dot{x}_i`` tracks the sensitivity to one of ``N``
independent inputs, and the ``\epsilon_i`` are formal infinitesimals with
``\epsilon_i\epsilon_j = 0``. Evaluating a function ``f`` on such a number yields

```math
f\!\left(x + \textstyle\sum_i \dot{x}_i \epsilon_i\right)
= f(x) + \sum_{i=1}^{N} \frac{\partial f}{\partial x_i}\,\dot{x}_i\,\epsilon_i,
```

so the value **and all ``N`` partials** come out of a single pass. Thus
"``N`` partials" means differentiating with respect to ``N`` inputs at once: a
turbomachinery cycle with ``N = 12`` design variables seeds 12 partials, and the
full gradient ``\nabla f \in \RR^{12}`` falls out of one evaluation.

### Why the analytic rules matter: flat in ``N``

Generic dual arithmetic carries the length-``N`` partials tuple through *every*
elementary operation in the NASA-9 polynomials, so the per-evaluation cost grows
roughly linearly in ``N`` (for ``dh/dT``: ``\approx 11`` ns at ``N=1`` rising to
``\approx 31`` ns at ``N=12``). The extension instead computes the single scalar
derivative once and scales the whole partials tuple,

```math
h\!\left(\text{gas},\, x + \dot{x}\right) = h(\text{gas}, x) + c_p(\text{gas}, x)\,\dot{x},
```

— one scalar times an ``N``-vector — so the cost is essentially **flat in
``N``** (``\approx 11`` ns at ``N = 1, 8, 12`` alike). See
[Performance](@ref performance) for the measured comparison.

## Inversions by the implicit function theorem

The inversions (`T_from_h`, and the polytropic engine behind `compress`/`expand`)
are bounded Newton solves. Differentiating the *iterates* would be expensive and
would leak solver tolerance into the gradient. Instead we use the **implicit
function theorem**: for a root ``T`` of ``F(T, p) = 0`` with parameters ``p``,

```math
\frac{dT}{dp} = -\left(\frac{\partial F}{\partial T}\right)^{-1}\frac{\partial F}{\partial p},
```

evaluated at the converged ``T^\star``. We solve for ``T^\star`` in plain
`Float64`, then attach the exact partials — the Newton loop never runs on dual
numbers.

For the **enthalpy inversion** ``F = h(\text{gas}, T) - h_\text{spec}``, we have
``\partial F/\partial T = c_p``, so

```math
dT = \frac{dh_\text{spec}}{c_p(\text{gas}, T^\star)}.
```

For the **polytropic engine** behind `compress`/`expand`, with pressure ratio
``\Pi`` and polytropic efficiency ``\eta_p``, the residual is
``G = s^0(T_2) - \big[s^0(T_1) + R\ln(\Pi)/\eta_p\big]``, giving

```math
\frac{c_p(T_2)}{T_2}\,dT_2 = \frac{c_p(T_1)}{T_1}\,dT_1 + \frac{R}{\eta_p \Pi}\,d\Pi .
```

This is about **2× faster** than differentiating the loop and, more importantly,
*exact*.

### Differentiating through composition

When the composition itself depends on a parameter — as in
`products(vit, FAR::Dual)`, where the burned-gas mole fractions are a function of
the fuel–air ratio — the resulting `FrozenGas` *carries* the tangent in its own
coefficients. The IFT rule then adds a "composition moves" term: the partials of
``h(\text{gas}, T^\star)`` with respect to the gas, obtained from a single forward
evaluation (the properties are linear in the coefficients), still without
touching the Newton loop.

## Worked example

Every block below runs when the docs are built.

```@example deriv
using IdealGasThermo, ForwardDiff

air = FrozenGas(DryAir)
T   = 800.0

# dh/dT recovers cp exactly
(dhdT = ForwardDiff.derivative(t -> h(air, t), T), cp = c_p(air, T))
```

```@example deriv
# dT/dh through the Newton inversion equals 1/cp (the IFT rule, exact)
hT = h(air, T)
(dTdh = ForwardDiff.derivative(hh -> T_from_h(air, hh), hT), inv_cp = 1 / c_p(air, T))
```

```@example deriv
# gradient of compressor-exit temperature w.r.t. (T₁, Π): N = 2 partials
f(x) = compress(GasState(air, x[1], 101325.0), x[2]; ηp = 0.9).T
ForwardDiff.gradient(f, [288.15, 12.0])
```

```@example deriv
# differentiate burned-gas enthalpy through combustion w.r.t. FAR —
# the composition depends on FAR, handled by the Dual-carrying-gas rule
vit = Vitiator("CH4", DryAir)
ForwardDiff.derivative(far -> h(products(vit, far), 1500.0), 0.03)
```
