# A Mathematical Theory of Aging in Humans and Dogs

10 November 2026

## What this paper proves and what it does not

This paper proves 19 theorems about aging. Together they say exactly what a cure for aging must do, how to recognise one in trial data, and why nothing available today is one. It does not contain a cure, and no theorem can.

The theorems rest on one measured fact: in adult humans and dogs the yearly risk of death rises by a fixed percentage each year. From that fact we prove:

- **Curing diseases is bounded.** Removing any cause of death adds at most a fixed number of years, set by how fast risk doubles (Theorems 2 to 4).
- **Dogs are not small, fast humans.** The gap between dog and human lifespans splits into two measurable parts, and one leading model of aging cannot turn a human into a dog by speed alone (Theorems 7 and 10).
- **Risk stops rising under one condition.** Damage must arrive at a pace that does not grow with age, and every kind of damage must be repairable (Theorem 12).
- **Escape velocity has an exact threshold.** A person stops aging in risk terms exactly when medicine lowers death rates each year as fast as their own risk grows: about 8 percent a year for humans, about 27 percent for dogs (Theorem 16).
- **One unrepaired kind of damage sets the limit.** Repairing everything else leaves risk doubling at a slower but nonzero pace (Theorem 18).

Each theorem comes with a proof and a computer check. The final sections turn the theorems into predictions that an experiment could confirm or kill. That step, not more mathematics, is what would make this science.

## Setup

Everything below follows from one assumption about how the risk of death changes with age, plus exact definitions of what "solving aging" could mean.

**Symbols**

| Symbol | Meaning |
| --- | --- |
| x | age, in years |
| t | calendar time, in years |
| μ(x) | hazard: the death rate per year among those alive at age x |
| S(x) | survival: the chance of living from birth to age x |
| e(x) | remaining life expectancy at age x |
| a | level: the aging part of the hazard at age 0 |
| b | rate of aging: the yearly growth of the hazard |
| c | background hazard that does not change with age |
| z(x) | the aging part of the hazard at age x, equal to a·e^(bx) |

**Assumption G (the Gompertz-Makeham law).**

```latex
\mu(x) = c + a\,e^{bx}, \qquad a > 0,\; b > 0,\; c \ge 0
```

The yearly risk of death is a constant background risk plus a part that grows by a fixed percentage every year.

This is a measured regularity, not a law of nature. It fits US adults aged about 30 to 95 and dogs aged about 2 to 17. It fails for infants and young adults, and it is disputed above age 105.

| Fitted from life tables | b per year | Risk doubles every |
| --- | --- | --- |
| Humans (US, 2023) | 0.083 to 0.102 | 6.8 to 8.3 years |
| Dogs (UK and US) | 0.31 to 0.34 | 2.0 to 2.3 years |

**Survival and life expectancy.**

```latex
S(x) = \exp\!\Big(-\int_0^x \mu(u)\,du\Big), \qquad e(x) = \int_0^\infty \exp\!\Big(-\int_0^s \mu(x+u)\,du\Big)\,ds
```

Survival falls with the total risk accumulated so far. Life expectancy at age x is the area under the survival curve from that age on.

**Definitions.** Four things "solving aging" could mean, weakest first.

1. **Slowed aging.** The rate b becomes smaller. Risk still rises with age, more slowly.
2. **Stopped aging.** The hazard no longer rises with age at any adult age.
3. **Escape velocity.** A person's own remaining life expectancy, recomputed each year from that year's death rates, never falls.
4. **Strong escape velocity.** At every fixed age, remaining life expectancy rises by at least one year per calendar year.

In symbols, with e(x, t) the remaining life expectancy at age x under the death rates of year t:

```latex
\text{(3)}\;\; \frac{d}{dt}\,e(x_0 + t,\; t) \ge 0 \qquad\qquad \text{(4)}\;\; \frac{\partial e}{\partial t}(x, t) \ge 1
```

Definition 3 follows one person as they age. Definition 4 compares people of the same age in successive years.

## Theorems 1 to 6: the mortality curve

Lowering the risk of death, by any amount short of 99 percent, cannot add more than a few doubling times of life. These six theorems make that exact.

### Lemma A (aging one year costs less than one year)

```latex
\frac{de}{dx} = \mu(x)\,e(x) - 1, \qquad\text{and if } \mu \text{ is increasing:}\quad -1 < \frac{de}{dx} < 0
```

Each year you age, your remaining life expectancy drops, but by less than a full year.

**Proof.** Write e(x) as the integral of S(y) from x to infinity, divided by S(x). Differentiating gives −1 + μ(x)e(x). If μ is increasing, survival from age x falls faster than exp(−μ(x)s), so e(x) is less than 1/μ(x), and μ(x)e(x) lies strictly between 0 and 1. ∎

### Theorem 1 (life expectancy in closed form)

Under Assumption G with c = 0, let M = z(x₀)/b be the aging hazard at age x₀ divided by the rate of aging. Then

```latex
e(x_0) = \frac{1}{b}\,e^{M}E_1(M), \qquad E_1(M) = \int_M^\infty \frac{e^{-u}}{u}\,du, \qquad e(x_0) \approx \frac{\ln(1/M) - \gamma}{b}\ \text{ for small } M
```

Remaining life depends on one number: how high your risk already is compared with how fast it is growing. Here γ = 0.5772 is Euler's constant.

**Proof.** Survival from x₀ is exp(−M(e^(bs) − 1)). Substituting u = M·e^(bs), so that ds = du/(bu), turns the integral for e(x₀) into (1/b)·e^M times the integral of e^(−u)/u from M to infinity. The approximation uses E₁(M) = −γ − ln M + O(M). ∎

**Numbers.** With the fitted values for US women (a = 3.03×10⁻⁵, b = 0.091) the formula gives 81.7 years at birth; the life table says 81.1.

### Theorem 2 (a risk cut is an age shift)

Under Assumption G, multiplying the aging hazard by a factor r between 0 and 1, from age x₀ onward, has exactly the same effect as making the person Δ years younger, where

```latex
\Delta = \frac{-\ln r}{b}, \qquad\text{and}\qquad 0 < \text{years gained} < \Delta
```

A treatment that cuts risk is worth a fixed number of years of youth, and never more.

**Proof.** r·a·e^(bx) = a·e^(b(x − Δ)), so the treated person's future hazard is that of an untreated person Δ years younger. The gain is therefore e(x₀ − Δ) − e(x₀), which is the integral of −de/dx over an interval of length Δ. By Lemma A that integrand lies strictly between 0 and 1. ∎

**Numbers.** A 20 percent cut (r = 0.8) from age 50 gains 2.24 years, against a bound of 2.45. For a 7-year-old dog the bound is 0.70 years.

### Theorem 3 (curing a disease is bounded)

Suppose total hazard obeys Assumption G with c = 0, causes of death act independently, and one cause never accounts for more than a share s̅ of the hazard at any age. Removing that cause entirely adds

```latex
\text{years gained} < \frac{-\ln(1 - \bar{s})}{b}
```

The gain from curing a disease is capped by its largest share of deaths, however deadly the disease.

**Proof.** After removal the hazard is (1 − s(x))μ(x), which is at least (1 − s̅)μ(x) at every age. More hazard means less survival at every age, so life expectancy is no greater than under a uniform cut with r = 1 − s̅. Apply Theorem 2. ∎

**Numbers (US, 2023).** Heart disease peaks at 26.7 percent of deaths (ages 85 and over): the bound is 3.1 to 3.7 years, and the full calculation gives 2.4 to 3.0. Cancer peaks at 23.1 percent (ages 45 to 64): the bound is 2.6 to 3.2 years, and the full calculation gives 2.3 to 2.7.

### Theorem 4 (halving all mortality buys less than one doubling time)

```latex
\text{years gained by halving every death rate} < \frac{\ln 2}{b}, \qquad \text{years gained} \ge Y \;\Rightarrow\; r < e^{-bY}
```

Halving every cause of death at every age adds under 8 years for a human and under 2.2 for a dog. Adding Y years by lowering risk alone needs a cut deeper than e^(−bY).

**Proof.** Set r = 1/2 in Theorem 2. For the second statement, the gain is below −ln(r)/b, so a gain of Y forces −ln(r)/b > Y. ∎

**Numbers.** Computed gains are 7.60 years (bound 7.62) for humans and 2.07 (bound 2.17) for dogs. Adding 50 human years needs r below 1.1 percent: a cut of about 99 percent in aging-related risk.

### Theorem 5 (the rate of aging matters far more than the level)

Under Assumption G with c = 0, hold the hazard at age x₀ fixed and divide the rate of aging by k > 1 from then on. Then

```latex
e_k(x_0) = \frac{k}{b}\,e^{kM}E_1(kM) \approx k\cdot\frac{\ln(1/M) - \ln k - \gamma}{b}
```

Slowing the rate of aging by a factor k multiplies remaining life by nearly k. Lowering the level by the same factor adds only ln(k)/b years.

**Proof.** Apply Theorem 1 with b replaced by b/k, so that M becomes kM. The comparison with the level follows from Theorem 2. ∎

**Numbers.** For a 30-year-old woman, halving the level adds about 7.5 years. Halving the rate adds about 37.

### Theorem 6 (how the record age grows)

Among N independent lives under Assumption G with c = 0, the median age of the oldest is

```latex
x_{\max} = \frac{1}{b}\ln\!\Big[1 + \frac{b}{a}\,\ln\frac{N}{\ln 2}\Big]
```

The record grows only with the logarithm of the logarithm of the population. If instead the hazard is constant at λ beyond some age, every tenfold increase in N adds ln(10)/λ years.

**Proof.** The oldest is below x with probability (1 − S(x))^N, which is one half when S(x) = ln(2)/N to first order. Solve for x using the Gompertz survival curve; for the constant-hazard tail, survival falls as e^(−λx). ∎

**Numbers.** For a billion women the formula gives 121.5 years; each tenfold increase in population adds about 1.5 years. If risk levels off at 0.645 per year after 105, as one Italian study found, each tenfold increase adds 3.6 years and there is no fixed ceiling. The verified record is 122 years 164 days.

## Theorems 7 to 9: dogs against humans

Dogs live about 7 times shorter lives than people, but age only about 3.5 times faster. These theorems split that gap into its parts and show what a result in dogs would mean for humans.

### Theorem 7 (why dogs live shorter: speed times level)

For two populations that obey Assumption G with c = 0 and small a/b, written h for humans and d for dogs,

```latex
\frac{e_h}{e_d} \approx \underbrace{\frac{b_d}{b_h}}_{\text{speed}} \times \underbrace{\frac{\ln(b_h/a_h) - \gamma}{\ln(b_d/a_d) - \gamma}}_{\text{level}}
```

The lifespan gap is the ratio of aging speeds times a second factor that measures how much lower the starting risk is in the longer-lived species.

**Proof.** Apply the approximation of Theorem 1 at birth to each population and divide. ∎

**Numbers.** With fitted values for US women and UK dogs the speed factor is 3.51 and the level factor is 2.07. Their product is 7.27; the exact ratio is 7.13. On a logarithmic scale, 63 percent of the gap is speed and 37 percent is level.

### Theorem 8 (dog years, and what a dog result means for humans)

Matching ages at which the same fraction of each population is still alive gives the map from dog age d to human age h:

```latex
h(d) = \frac{1}{b_h}\ln\!\Big[1 + \frac{a_d\,b_h}{a_h\,b_d}\big(e^{b_d d} - 1\big)\Big], \qquad \frac{dh}{dd} \to \frac{b_d}{b_h}
```

In adulthood every dog year is worth the same number of human years, the ratio of aging speeds, which is about 3.5.

**Corollary (translation rule).** If a treatment has the same hazard ratio r in both species, the years of youth it buys are in the ratio

```latex
\frac{\Delta_h}{\Delta_d} = \frac{b_d}{b_h}
```

Each dog year gained would correspond to about 3.5 human years. Nothing guarantees that the hazard ratio carries over; that is an empirical question.

**Proof.** Set the two Gompertz survival curves equal and solve for h. As d grows the 1 inside the logarithm becomes negligible, leaving a straight line of slope b\_d/b\_h. The corollary is Theorem 2 applied to each species. ∎

**Numbers.** The slope is 3.53 at dog age 16 and tends to 3.51. Real life tables, which include the risky first years, map a 1-year-old dog to a 33-year-old human and a 12-year-old dog to 78.

### Theorem 9 (why dog trials are faster)

A trial that compares treatment with placebo in equal groups needs a fixed number of deaths D to detect a hazard ratio r, and a number of participants N that depends on how fast deaths arrive:

```latex
D = \frac{4\,(z_{1-\alpha/2} + z_{1-\beta})^2}{(\ln r)^2}, \qquad N = \frac{D}{\bar{p}}, \qquad \bar{p} = 1 - \tfrac{1}{2}\big(e^{-H} + e^{-rH}\big), \qquad H = cT + \frac{z(x_0)}{b}\big(e^{bT} - 1\big)
```

The deaths needed depend only on the size of the effect. The people or dogs needed depend on the accumulated risk H over T years of follow-up, which grows like e^(bT).

**Proof.** The first formula is Schoenfeld's (1981) and is taken as given. The expression for H is the integral of the hazard in Assumption G from x₀ to x₀ + T; the chance of dying in a group is one minus the exponential of minus its accumulated hazard. ∎

**Consequence.** When background risk is small and both species enrol at matching ages (the same z(x₀)/b), H depends on T only through bT, so a species that ages k times faster reaches the same accumulated risk in one k-th of the time. Dogs compress a human trial by about 3.5.

**Numbers.** Detecting a 20 percent cut in risk at the usual standards needs 631 deaths. That is about 4,000 large dogs aged 7 followed for 3 years, or about 16,700 people aged 65.

## Theorems 10 to 14: mechanisms

These theorems concern models of what goes wrong inside the body. They are exact statements about the models. Whether a model describes a real body is a separate, experimental question.

### The worn-out-cell model

Damaged "senescent" cells X are produced at a rate that grows with age t, cleared by a system that saturates, and jostled by noise. Death occurs when X first reaches a threshold X\_c (Karin and others, 2019).

```latex
dX = \Big[\eta\,t - \frac{\beta X}{\kappa + X}\Big]dt + \sqrt{2\epsilon}\;dW
```

Here η is how fast production climbs, β is the most the body can clear, κ is where clearing starts to saturate, and ε is the size of the random swings. The cell load settles in days to weeks while production changes over years, so at each age the hazard is the hazard for today's production rate p = ηt held fixed. Call it H(p).

### Theorem 10 (speed alone cannot turn a human into a dog)

In the worn-out-cell model, multiplying the production slope η by s rescales age:

```latex
h(t;\, s\eta) = h(st;\, \eta)
```

If the hazard is close to a·e^(bt) over adult ages, the faster species has rate of aging s·b and the same level a. Its median lifespan L then satisfies

```latex
(sb)\,L_s \approx b\,L + \ln s
```

Speeding up production multiplies the rate of aging but leaves the level untouched, so rate times lifespan can only go up.

**Proof.** The hazard is H(ηt), so replacing η by sη gives H(sηt) = h(st; η). Under a·e^(bt) this is a·e^(sbt): same a, rate sb. The Gompertz median is ln(1 + b·ln2/a)/b, so rate times median is ln(b·ln2/a) to leading order, and replacing b by sb adds ln s. ∎

**Numbers.** Rate times median lifespan is about 7.6 for humans and about 3.8 for dogs. A human sped up 3.5 times through production alone would be a "dog" with a median lifespan near 28 years, not 12. Fitted levels differ by a factor of about 160. So in this model dogs must also clear worn-out cells more weakly, have a lower lethal threshold, or suffer larger random swings.

### Theorem 11 (what sets the rate of aging in this model)

When the lethal threshold is far above the typical cell load, the hazard grows at the rate

```latex
\frac{d\ln h}{dt} \approx \frac{\eta\,(X_c - X^*)}{\epsilon} - \frac{\eta}{U'(X_c)} - \frac{\eta}{\beta - p}, \qquad X^* = \frac{\kappa p}{\beta - p}
```

The rate of aging is the climb in production, times the distance from the usual load X\* to the lethal level, divided by the noise. The two smaller terms slow the rise late in life. U is the potential whose slope is clearance minus production.

**Proof sketch.** For a high barrier the escape rate is a prefactor times exp(−ΔU/ε), where ΔU is the height of the barrier between X\* and X\_c. Differentiating ΔU with respect to p gives −(X\_c − X\*), because the term from the moving minimum vanishes. The other two terms are the derivative of the logarithm of the prefactor. This is an asymptotic result, not an identity. ∎

**Numbers.** With the published human settings the formula gives 0.130, 0.116 and 0.081 per year at ages 40, 60 and 80. Direct computation gives 0.128, 0.114 and 0.090.

### Theorem 12 (when risk stops rising)

Consider any model in which the body moves among finitely many living states at rates that do not depend on age, every living state can reach every other, and death is possible. Then the hazard among survivors converges to a constant:

```latex
\lim_{t\to\infty} h(t) = \lambda, \qquad -\lambda = \text{largest eigenvalue of the matrix } Q \text{ of rates among living states}
```

If nothing in the damage process changes with age and all damage can be undone, risk levels off and stays flat.

**Proof.** The distribution over living states is p(t) = p(0)·e^(Qt). Q has non-negative off-diagonal entries and is irreducible, so by the Perron-Frobenius theorem p(t) approaches a multiple of e^(−λt)·v with v positive and vQ = −λv. The hazard is p(t)·d divided by the total of p(t), where d is the vector of death rates and Q·1 = −d. Hence the limit is v·d divided by the total of v, which equals λ. ∎

**Corollary.** A risk that keeps rising without limit needs one of two things: rates that change with age, or damage that can pile up without bound because it is never repaired. Those are the only two targets for stopping aging.

**Caution.** The theorem does not say how long levelling off takes. It can take longer than a lifetime.

### Theorem 13 (a uniformly dangerous world does not select for faster aging)

Let l(x) be survival to age x, m(x) fertility, and r the population growth rate, which solves

```latex
\int_0^\infty e^{-rx}\,l(x)\,m(x)\,dx = 1
```

Adding the same extra death rate at every age changes r but leaves unchanged the strength of selection against dying at each age. That strength is zero after the last age of reproduction.

**Proof.** An extra constant hazard δ multiplies l(x) by e^(−δx). The equation is then solved by r − δ, so e^(−rx)·l(x) is unchanged at every age. Hamilton's measure of selection against death at age x is the integral of e^(−ry)·l(y)·m(y) over ages y beyond x, divided by generation time. Both are built from that unchanged product, and the integral is zero once m is zero. ∎

**Meaning.** Evolution lets aging happen because it cannot see what occurs after reproduction ends. Whether a dangerous environment speeds aging depends on how population size is held in check, so that popular claim is not automatic.

### Theorem 14 (cancer risk peaks in medium-large dogs)

Suppose lifetime cancer risk is proportional to body weight W times lifespan T to the power n, and lifespan falls with weight as T = T₀ − kW. Then risk is greatest at

```latex
W^* = \frac{T_0}{(n+1)\,k}
```

More cells raise the risk; a shorter life cuts it back; the peak sits where the two balance.

**Proof.** Differentiate W·(T₀ − kW)^n and set the result to zero: T₀ − kW = nkW. ∎

**Numbers.** With the values fitted to 85 breeds (n = 4, T₀ = 13.3 years, k = 0.08 years per kg; Nunney 2024), the peak is at 33.25 kg.

## Theorems 15 to 19: escape velocity and repair

Medicine outruns aging exactly when it lowers death rates each year as fast as a person's own risk grows. For humans that is about 8 percent a year, against an observed 0.5 to 1 percent.

Throughout this section medicine improves steadily: the aging part of every death rate falls by a constant fraction ρ per year.

```latex
\mu(x, t) = c + a_0\,e^{-\rho t}\,e^{bx}
```

In calendar year t, a person aged x faces the usual law, with the level lowered by t years of progress.

### Theorem 15 (progress is an age shift)

```latex
\frac{\partial e}{\partial t} = -\frac{\rho}{b}\,\frac{\partial e}{\partial x} = \frac{\rho}{b}\,\big(1 - \mu e\big)
```

Lowering aging-related death rates by a fraction ρ each year makes everyone ρ/b years younger each year.

**Proof.** Life expectancy depends on x and t only through the combination z = a₀·e^(bx − ρt). By the chain rule its derivative in t is −ρ/b times its derivative in x. The second equality is Lemma A. ∎

### Theorem 16 (the escape velocity threshold)

Follow a person aged x₀ at time 0. Their own remaining life expectancy never falls if and only if

```latex
\rho \ge b
```

At every fixed age, remaining life expectancy rises by at least one year per year if and only if

```latex
\rho \ge \frac{b}{1 - \mu e}
```

You stop aging in risk terms when medicine improves as fast as you age. The stricter version needs faster progress, and more at older ages.

**Proof.** Along a life, age and calendar time advance together, so the rate of change is the sum of the two partial derivatives. By Theorem 15 that sum is (1 − ρ/b)(μe − 1). Lemma A makes the second factor negative, so the product is non-negative exactly when ρ ≥ b. The second statement sets the right-hand side of Theorem 15 to at least 1. ∎

**Remark.** Along that life the hazard is c + z₀·e^((b − ρ)t). At ρ = b it is constant: the person is frozen at their current risk, not made immortal.

**Numbers.** A continuous rate ρ = b means death rates fall 8.1 to 9.5 percent a year for humans and 27 percent for dogs. US life expectancy at 65 rose from 13.9 years (1950) to 19.7 (2024), which corresponds to 0.5 to 1.0 percent a year: 8 to 17 times too slow. Remaining years for a 50-year-old at 0, 50, 90, 100 and 110 percent of the threshold: 33, 52, 123, 272 and about 1,320.

### Theorem 17 (periodic repair)

Suppose damage D builds up between treatments and each treatment, given every T years, removes a fraction f of it.

If damage builds at a steady rate r, then for any f above zero it settles into a fixed cycle:

```latex
D_{\text{after}} = rT\,\frac{1-f}{f}, \qquad D_{\text{before}} = \frac{rT}{f}
```

If damage also feeds on itself at rate g, it stays bounded if and only if

```latex
f > 1 - e^{-gT}
```

Steady damage can be held in check by any regular partial repair. Self-feeding damage must be removed faster than it regrows.

**Proof.** One cycle maps D to (1 − f)(D + rT). This is a contraction with factor 1 − f, and its fixed point is the stated value. With self-feeding growth one cycle maps D to (1 − f)(D·e^(gT) + (r/g)(e^(gT) − 1)), which contracts exactly when (1 − f)·e^(gT) is below 1. ∎

**Numbers.** Removing half the damage every 10 years holds a body between 10 and 20 years' worth of damage indefinitely.

### Theorem 18 (the weakest link sets the limit)

Let the logarithm of the aging hazard be a sum of K kinds of damage, growing at rates r₁ to r\_K that add up to b. Suppose repair holds every kind bounded except kind j. Then along a life

```latex
a\,e^{r_j t} \;\le\; z(t) \;\le\; C\,a\,e^{r_j t}
```

for a constant C. The rate of aging falls from b to r\_j, and risk doubles every ln(2)/r\_j years. Aging stops only if every kind of damage is held bounded.

**Proof.** The repaired kinds contribute a bounded amount to the logarithm of the hazard, between 0 and ln C. The unrepaired kind contributes r\_j·t. ∎

**Numbers.** With seven equal kinds, repairing six perfectly and never the seventh leaves risk doubling every 57 years in place of every 8. In an illustrative model that gives an expected age at death near 380, against about 1,800 with all seven repaired.

### Theorem 19 (population growth is set by births)

The long-run growth rate of a population depends on death rates only up to the last age of reproduction. Changing mortality after that age leaves the growth rate exactly unchanged. A population that has stopped growing has size

```latex
N = (\text{births per year}) \times e(0)
```

Ending aging makes a population larger without making it grow faster.

**Proof.** The growth rate solves the equation of Theorem 13, whose integrand is zero wherever fertility is zero. With zero growth, the number alive at each age is births per year times survival to that age; summing over ages gives births times life expectancy. ∎

**Numbers.** At 1.5 children per woman a population shrinks about 1.1 percent a year in the long run, with or without aging.

## Predictions an experiment can confirm or kill

The theorems become science only where they stick their neck out. Each prediction below names a measurement and the result that would prove it wrong.

### Prediction 1: dogs clear worn-out cells differently, not just faster

- **From:** Theorem 10.
- **Claim:** If the worn-out-cell model describes both species, a young adult dog must have weaker clearance of senescent cells, a lower lethal load, or larger random swings than a human sped up 3.5 times.
- **Test:** Measure how fast senescent cells are cleared in young and old dogs, as has been done in mice, and compare with human tissue.
- **Killed if:** Dog cell dynamics are simply human dynamics run 3.5 times faster. The model would then fail to explain why young dogs carry about 160 times the fitted level of risk.

### Prediction 2: a treatment's fingerprint shows whether it touched aging

- **From:** Theorems 2, 4 and 5.
- **Claim:** A treatment that only lowers risk shifts the curve of log death rate against age sideways, keeping its slope. A treatment that slows aging changes the slope. A level-only treatment can add at most one doubling time unless it cuts risk by more than half: 2.2 years in dogs.
- **Test:** In the dog trials now running (1,300 and 580 dogs), plot log death rate by age for treated and placebo groups and compare slopes.
- **Killed if:** A trial shows a gain beyond one doubling time with parallel curves and a risk cut under half. That would contradict Assumption G for treated animals.

### Prediction 3: one dog year gained is worth about 3.5 human years

- **From:** Theorem 8.
- **Claim:** If a drug's hazard ratio is the same across species, the years it adds scale with the inverse of the rate of aging.
- **Test:** Compare hazard ratios for one drug in mice, dogs and people. Rapamycin is the first candidate: it has mouse results and a dog trial under way.
- **Killed if:** Hazard ratios differ systematically between species. The scaling rule would then have no practical use.

### Prediction 4: a valid aging marker must predict deaths in finished trials

- **From:** Theorem 2, applied to a marker.
- **Claim:** If a blood or DNA marker truly stands in for aging, then a treatment that shifts the marker by d standard deviations must change the risk of death by the factor

```latex
r = h^{\,d}
```

Here h is the hazard ratio per standard deviation of the marker in untreated people. The marker's movement must account for the change in deaths.

- **Test:** Measure markers in stored blood from finished trials that counted deaths, such as the large blood-pressure and semaglutide trials.
- **Killed if:** The change in deaths is not predicted by the change in the marker. That marker cannot replace lifespan in trials.

### Prediction 5: the record age tells us whether risk levels off

- **From:** Theorems 6 and 12.
- **Claim:** If risk keeps doubling, each tenfold increase in the number of people reaching extreme age raises the record by about 1.5 years. If risk levels off near 0.65 per year, it raises it by about 3.6 years.
- **Test:** Track the oldest verified ages as the international database of people aged 105 and over grows.
- **Killed if:** Neither pattern holds, which would point to errors in age records.

**What none of these would do.** Confirming all five would validate the theory as a tool. It would not produce a treatment. Predictions 2 and 4 would, however, tell the field quickly whether a given treatment is touching aging itself.

## Checks, what is new, limits and references

Every theorem passed a computer check, and most of them restate or sharpen results that are already known. The honest claim for this paper is a unified, verified framework with five testable predictions, not a discovery.

### Computer checks

One script, `verify_theorems.py`, tests each theorem numerically or symbolically and stops at the first failure. It runs in under a minute and is in the `theory/` folder of this repository.

| Theorem | Check | Result |
| --- | --- | --- |
| 1 | Closed form against direct integration, 4 cases | agree to 8 digits |
| 2 | Gain below the bound, and the age-shift identity, 20 cases with and without background risk | holds in all; identity exact |
| 3 | 200 random diseases with random share curves | gain never above 94% of the bound |
| 4 | Halving all mortality | 7.60 < 7.62 years (human), 2.07 < 2.17 (dog) |
| 5 | Halving the rate against halving the level | +37.2 against +7.5 years |
| 6 | Record-age formula against a root-finder | agree to 6 digits |
| 7 | Speed times level against the exact ratio | 7.27 against 7.13 |
| 8 | Survival match and limiting slope | exact; slope 3.511 |
| 10 | Rescaling identity; rate times lifespan shifts by ln s | identity to 8 digits; 8.90 predicted and exact |
| 11 | Formula against direct computation | within 2% at ages 40 to 60 (asymptotic) |
| 12 | Three-state repairable chain | hazard converges to the eigenvalue, 0.0844 |
| 13, 14 | Symbolic algebra | confirmed |
| 15, 16 | Both identities; life expectancy along a life at ρ = b | agree to 9 digits; constant at 31.04 years |
| 17 | Iterating 200 to 400 treatment cycles | fixed point exact; threshold 0.259 confirmed |
| 18 | Long-run slope with one unrepaired kind of seven | 0.012145 against 0.012143 predicted |
| 9, 19 | Recomputed independently by a second agent | confirmed |

### What is classical and what may be new

**Classical, restated here.** The mortality law (Gompertz 1825, Makeham 1860). Lemma A and the closed form of Theorem 1. The small gain from curing a disease, known as the Taeuber paradox (Keyfitz 1977). The record-age formula. Schoenfeld's trial formula. Theorem 12, which is standard Markov-chain theory and the basis of reliability models of aging (Gavrilov and Gavrilova). Theorem 13 (Hamilton 1966; Moorad and others 2019). Theorem 14 (Nunney 2024). Theorem 19 (Lotka). The worn-out-cell model (Karin and others 2019). The idea of escape velocity (de Grey 2004).

**Possibly new, or newly assembled.** The strict bounds of Theorems 2 to 4 in this form. The speed-and-level split of the dog and human gap (Theorem 7). The argument that speed alone cannot produce a dog (Theorem 10). The escape threshold with a background term (Theorem 16). The weakest-link bound (Theorem 18). We have not done a full search of the literature, so none of these should be claimed as original until an expert has checked.

### Limits

- Everything rests on Assumption G, which is a fit to data from mid-adulthood to old age.
- Theorems 10 and 11 are about one model of one mechanism. Its settings cannot be pinned down from death records alone.
- Theorems 15 to 18 assume steady, uniform progress and idealised repair.
- No theorem shows that any treatment works. The mathematics says what a cure must do and how to recognise one.

### What would make this publishable

1. An expert in the mathematics of aging reviews the proofs and the claims of novelty.
2. One prediction is tested. Predictions 2 and 4 need only data that exist or are being collected now.
3. The result, whichever way it falls, is written up with the theorem it tests.

### References

- [US Social Security period life table, 2023](https://www.ssa.gov/oact/STATS/table4c6.html)
- [CDC, Deaths: Leading Causes for 2023](https://cdc.gov/nchs/data/nvsr/nvsr74/nvsr74-10.pdf)
- [Teng 2022, UK dog life tables](https://pmc.ncbi.nlm.nih.gov/articles/PMC9050668)
- [Montoya 2023, US dog life tables](https://www.frontiersin.org/journals/veterinary-science/articles/10.3389/fvets.2023.1082102/pdf)
- [Keyfitz 1977, what difference would it make if cancer were eradicated](https://gwern.net/doc/longevity/1977-keyfitz.pdf)
- [Barbi 2018, death rates after age 105](https://pmc.ncbi.nlm.nih.gov/articles/PMC6457902)
- [Gavrilov and Gavrilova 2001, reliability theory of aging](https://researchgate.net/profile/Natalia-Gavrilova-5/publication/11613344_The_Reliability_Theory_of_Aging_and_Longevity/links/653b472f24bbe32d9a7216d9/The-Reliability-Theory-of-Aging-and-Longevity.pdf)
- [Karin 2019, senescent cell turnover and the Gompertz law](https://pmc.ncbi.nlm.nih.gov/articles/PMC6889273)
- [Moorad, Promislow and Silvertown 2019, evolution of senescence](https://www.zoology.ubc.ca/~otto/veg/Readings/Moorad2019.pdf)
- [Nunney 2024, body size and cancer in dog breeds](https://escholarship.org/content/qt8w0923ff/qt8w0923ff.pdf?t=s8hg2k)
- [Kraus, Pavard and Promislow 2013, why large dogs die young](https://jhanley.biostat.mcgill.ca/bios601/Applications/LongevityOfDogs.pdf)
- [Oeppen and Vaupel 2002, record life expectancy](https://u.demog.berkeley.edu/~jrw/Biblio/Eprints/%20M-O/oeppen.vaupel.2002.pdf)
- [Olshansky 2024, limits to life expectancy](https://www.nature.com/articles/s43587-024-00702-3)
