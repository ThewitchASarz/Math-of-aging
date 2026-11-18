# Experiment 017: result

**The hidden-differences slope test of experiment 012 holds up in a 4-year trial and in a 10-year trial with dropout and a mix of sizes: a plain 50 percent cut in risk is flagged as a change in the rate of aging in about 4 percent of trials, and no power is lost. It does not hold up when the treatment helps some dogs more than others. If frail dogs are helped more, the new test flags 30 percent of 10-year trials and the plain test 78 percent, always as faster aging. If sturdy dogs are helped more, the new test flags 88 percent of trials and the plain test 58 percent as slowed aging, although no dog's rate of aging changed.** Two predictions were confirmed and two were half confirmed and half refuted.

This is a calculation under stated assumptions. It says how two analyses behave on simulated dogs, not that any treatment works. Nobody has measured how widely dogs differ in hidden risk, or whether any treatment helps frail and sturdy dogs differently; the spread of 0.3 and the 70 and 30 percent cuts are assumptions.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | 4-year trial, 50% risk cut: plain test flags 6.5 to 15%, new test 3 to 8%; new test with no effect 3 to 8% | Plain 9.1% (check: 8.7%); new 4.1% (check: 4.3%); no effect 5.2% (check: 4.3%) | Confirmed |
| P2 | 4-year trial, aging slowed 20%: new test flags 15 to 32%; with no hidden differences the new test is no more than 5 points below the plain test | 25.8% (check: 25.8%); with no hidden differences, new 29.4% against plain 28.4% | Confirmed |
| P3 | 10 years with sizes and dropout, 50% risk cut: plain test flags 12 to 25%, new test 3 to 8%; new test with no effect 3 to 8%; new test recognises a 20% slowing in 60 to 80% | Plain 8.7% (check: 12.3%); new 3.6% (check: 4.6%); no effect 4.1%; slowing 80.5% (check: 77.2%) | **Half and half.** The new test behaved as predicted. The plain test was disturbed less than predicted in our run, and the power figure is on the upper line |
| P4 | Frail dogs helped more, 10 years: plain test flags 25 to 60%, new test more than 10%; fewer than 3% of trials flagged as slowed aging | Plain 77.6% (check: 77.7%); new 30.5% (check: 30.5%); none flagged as slowed aging | **Half and half.** Both tests are defeated and the error points away from "slowed aging", as predicted. The plain test is defeated far more often than predicted |

## Main table

Share of 800 simulated trials in which each test reports a change in the rate of aging. Every trial has 1,300 dogs, half treated from enrolment, and dogs with widely spread hidden risk (spread 0.3, ASSUMPTION) unless stated. Each figure is uncertain by about 1 to 1.5 percentage points.

| Trial | True effect | Plain test | New test | Should be |
| --- | --- | --- | --- | --- |
| 4 years, entry ages 10 to 14 | none | 5.4% | 5.2% | about 5% |
| | 50% risk cut | 9.1% | 4.1% | about 5% |
| | aging slowed 20% | 23.1% | 25.8% | high |
| | aging slowed 20%, no hidden differences | 28.4% | 29.4% | high |
| 10 years, five sizes, 5% dropout a year | none | 6.4% | 4.1% | about 5% |
| | 50% risk cut | 8.7% | 3.6% | about 5% |
| | aging slowed 20% | 70.6% | 80.5% | high |
| 10 years, treatment helps frail dogs more (cuts of 70% and 30%) | no change in rate | **77.6%**, all as faster aging | **30.5%**, all as faster aging | about 5% |
| 10 years, treatment helps sturdy dogs more (cuts of 30% and 70%) | no change in rate | **58.0%**, all as slowed aging | **88.0%**, all as slowed aging | about 5% |
| 4 years, treatment helps frail dogs more | no change in rate | **24.0%**, all as faster aging | 5.1% | about 5% |

## What this means

1. **Three of the four limits on experiment 012's repair can be eased.** The extra setting works in a 4-year trial in which 60 to 70 percent of dogs die, and in a 10-year trial with five sizes of dog and about a quarter of dogs lost to follow-up. In both it brings false alarms for a plain risk cut to about 4 to 5 percent, and it recognises a true slowing at least as often as the plain test.
2. **A 4-year trial still has little power to recognise slowed aging**, with or without the extra setting: about 23 to 29 percent for a 20 percent slowing. That repeats experiments 001 and 016. The repair makes the test honest; it does not make a short trial long.
3. **Both tests assume the treatment helps every dog by the same fraction. When it does not, both can be badly wrong.** In each dog the benefit is a fixed cut in risk. But the dogs that die early are mostly the frail ones and the dogs that die late are mostly the sturdy ones. So the gap between the groups early in the trial shows the benefit to frail dogs, and the gap late in the trial shows the benefit to sturdy dogs. If the two benefits differ, the gap drifts over time, and a drifting gap is exactly what the slope test reads as a change in the rate of aging.
4. **The dangerous direction is a treatment that helps sturdy dogs more.** Then the gap between the groups widens with time, and both tests report slowed aging in most trials (58 and 88 percent) with an estimated slowing of 20 to 25 percent, when no dog aged more slowly. We said beforehand that such a result would go at the top of `KNOWN_GAPS.md`, and it has. This was the mirror case, reported and not scored.
5. **The new test is not always the safer one.** It is safer when frail dogs are helped more (30 against 78 percent, and 5 against 24 percent in a 4-year trial). It is worse when sturdy dogs are helped more (88 against 58 percent). Why it is worse there was not established; the check's suggestion is that the new test then estimates a small spread (about 0.1 where the truth is 0.3) and loses the caution the spread normally gives it.

For a real trial this suggests: a finding of "slowed aging" from a slope test, with or without a hidden-differences setting, also fits a treatment that simply helps the sturdiest dogs most. A trial would need some measured sign of frailty at entry, and a planned look at whether the benefit differs by it, to tell the two apart. Whether that look would work has not been tested here. This is a statement about the mathematics, not a comment on any trial's own analysis plan, which answers its own questions.

## Independent check

A separate agent derived the survival formula itself and agreed with ours. It wrote its own simulation and fitting code from the description before reading ours, and ran all ten cases with its own seeds (400 to 800 trials per case). Its figures are in the table above. Eighteen of twenty comparisons were within two standard errors; the two largest gaps were the plain test in the 10-year sizes-and-dropout trial with a risk cut (8.7 percent here, 12.3 in the check's 800 trials) and the plain test in the 4-year trial with no effect (5.4 here, 3.0 in the check's 400). The share of dogs dying agreed in every case to three digits.

It then read our code. It ran our functions and its own on the same simulated trials: the plain test's statistic agreed to better than 0.001 and every decision matched; where the new test's statistic differed, a brute-force search sided with our fit. It checked our gradient numerically (error about 1 part in 100 million), found no statistic below zero in 300 trials and no estimate stuck on a limit other than the spread, and confirmed the simulation (selection of frail dogs at entry, the frail and sturdy halves, dropout, alternating groups), the adding-up of results, and every input against its file in this repository.

It reproduced the sturdy-helped-more result (plain 55.7 percent, new 88.0 percent, all as slowed aging) and judged it a real feature of the model, not a coding fault. It also computed, with no chance involved, the ratio of death risk between treated and untreated dogs alive at each time, for dogs entering at age 9.5:

| Years in the trial | Frail helped more | Sturdy helped more |
| --- | --- | --- |
| 0.5 | 0.42 | 0.58 |
| 2 | 0.44 | 0.57 |
| 4 | 0.50 | 0.54 |
| 6 | 0.58 | 0.48 |
| 8 | 0.73 | 0.43 |

The ratio rises when frail dogs are helped more and falls when sturdy dogs are helped more. No single dog's ratio changes.

What the check changed:

- **The plain test's false-alarm rate in the 10-year sizes-and-dropout trial is stated as "about 9 to 12 percent"**, not 8.7. The gap between the two runs is about 2.3 standard errors, so part of our miss of P3's range may be chance. The verdict on P3 stays half and half.
- **A comment in the code was wrong** (it said the slowest cases ran first; the order was different). It only affects scheduling. Corrected.
- It noted that our plain fits of the larger model do not start from the smaller model's full answer, so "the larger model never fits worse" is not guaranteed by construction. It held in every trial tested.
- It added limit 5 below, and the reason given in point 5 above is its suggestion, not an established fact.

**Not checked:** spreads other than 0.3; whether the estimated spread is itself accurate; our share of dogs dropping out (about 23 to 25 percent) beyond a plausibility look. The checker's own fits of the new test stopped slightly early when the spread reached its lower limit, so its new-test figures in two cases are slight underestimates.

## How this compares with published work

That hidden differences between individuals bend a mortality curve, and that survivors of two trial groups can differ in their mix, is published: `papers/literature-check.md` credits the first to Beard (1959) and Vaupel, Manton and Stallard (1979), read there through Missov et al. (2015), and the second to Bartlett et al. (2020), who show that risk ratios for separate periods are biased by this selection even in a randomised trial. The drifting ratio found here is the same kind of selection effect, with a benefit that differs between dogs added. The literature check found no hidden-differences model fitted to dogs, no measured spread for dogs, and no published error rates for a slope test under any of these conditions in the searches run, so there is no published figure to compare these numbers with. Work on treatment benefits that differ between frail and sturdy individuals was not searched for; that is a gap in the search, not a claim that none exists.

## Limits

1. **ASSUMPTION: a spread of 0.3.** No measurement exists. The only fitted value in this repository is about 0.1 (experiment 004), from a table that mixes breeds. With a smaller spread every effect here would be smaller.
2. **ASSUMPTION: the 70 and 30 percent cuts, and two halves.** A real difference in benefit could be smaller, could vary smoothly, or could be absent. One size of difference was tried, in one design each way. The 4-year and sizes-and-dropout trials were not run with the sturdy-helped-more case.
3. **ASSUMPTIONS: entry ages (10 to 14 and 7 to 12), dropout of 5 percent a year unrelated to health, the US size mix, and the same spread inside every size class.** Experiment 014 showed dropout linked to health is a separate problem.
4. **The analysis has the right shape for the simulated dogs in Parts 1 and 2**: a doubling law, a gamma-shaped hidden factor, an effect constant from enrolment. Experiment 012 tried one other shape; this experiment tried none.
5. **In a 4-year trial the spread is poorly estimated.** The estimate sits at its lower limit in 16 to 28 percent of trials when the true spread is 0.3. The test still behaved, but the estimated spread itself should not be reported as a finding from a short trial.
6. **One trial size (1,300 dogs), one size of risk cut (50 percent), one size of slowing (20 percent).**
7. With 800 trials per case, differences of 2 to 3 points between cases can be chance.

## Inputs

- Level and slope of dog mortality: `dog_lifespan/results.json` (fit to the UK life table of Teng et al. 2022).
- Size-class rates of aging: `dog_lifespan/data/size_class_gompertz_calibration.csv`, scenario A; size mix from `dog_lifespan/data/montoya2023_size_e0.csv` (Montoya et al. 2023). Both read through experiment 007's code.
- Likelihoods: experiment 001 (`fingerprint.py`), experiment 007 (`breeds.py`), experiment 012 (`hidden_differences.py`).
- No outside source was opened in this run. Raz et al. 2026 and McMillan et al. 2024 were tried again for queue items 022 and 023 and could not be opened.

## Files

`harder.py` (about 5 minutes), `results.json`, `PREDICTION.md`.
