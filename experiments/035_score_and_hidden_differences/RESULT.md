# Experiment 035: result

**Adding a setting for left-over hidden differences to the score-adjusted slope test makes the unequal-benefit trap worse, not better. With sturdy dogs helped more and a frailty score of quality 0.7, the test with the extra setting reports slowed aging in about 59 percent of simulated trials, against about 40 percent without it, though no dog aged more slowly. With no score at all it is 87 percent against 52. Yet the same setting is still needed for the opposite mistake: without it, a plain 50 percent cut in risk is reported as faster aging in about 24 percent of trials with no score and about 11 percent with a score of quality 0.7. So no single version of the test is safe against both mistakes.** All four predictions were confirmed.

This is a calculation under stated assumptions. It says how analyses behave on simulated dogs. It does not say that any treatment works, that any treatment helps sturdy dogs more, or that any real frailty measurement reaches any of these qualities. The spread of hidden risk (0.3), the form of the scores and every size of benefit are assumptions.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | No effect: the test with the spread setting flags 3 to 8% of trials with no score, with the rank score at quality 0.7 and at 1.0 | 5.2%, 4.3% and 5.0% (check: 4.3%, 4.7% and 4.7%) | Confirmed |
| P2 | Sturdy helped more: the test with the spread setting reports slowed aging in 80 to 95% with no score, 48 to 70% at quality 0.7 and 24 to 42% at quality 0.9 | 87.2% (check: 88.0%); 59.3% (check: 59.7%); 27.5% (check: 24.9%) | Confirmed |
| P3 | Risk cut 50% in every dog: without the spread setting, faster aging is reported in 12 to 28% with no score and 8 to 22% at quality 0.7; with it, 3 to 9% are flagged in both | 24.2% and 11.5% (check: 24.5% and 11.0%); 5.0% and 5.3% (check: 5.3% and 7.1%) | Confirmed |
| P4 | Aging truly slowed 20%: the test with the spread setting recognises it in 55 to 74% at quality 0.7, no more than 8 points below the test without it, and in 58 to 80% with no score | 67.5% against 62.2% without (check: 67.9% against 62.2%); 69.7% (check: 68.5%) | Confirmed |

## The two analyses

Both are the slope test of experiment 034. It includes a frailty score measured at entry, lets the benefit differ by score in a straight line, and gives entry age its own setting. Each asks: did the rate at which risk rises change in treated dogs?

- **Without spread**: experiment 034's test. It assumes the score captures every difference between dogs.
- **With spread**: one more setting for hidden differences the score does not see. This is experiment 012's idea, reworked so that entry age keeps its own setting.

"No score" means the analysis has no score terms at all. The rank score and the risk score are those of experiment 034. Quality runs from 0 (useless) to 1 (perfect).

## Main tables

Each figure is the share of 600 simulated trials of 1,300 dogs followed for 10 years. Each is uncertain by about 1 to 2 percentage points, so differences of 2 to 4 points mean little.

**Treatment does nothing.** Share of trials flagged as a change in the rate of aging. Should be about 5 percent.

| Score | Without spread | With spread |
| --- | --- | --- |
| No score | 5.0% | 5.2% |
| Rank, quality 0.7 | 5.2% | 4.3% |
| Rank, quality 0.9 | 4.7% | 4.3% |
| Rank, quality 1.0 | 5.3% | 5.0% |
| Risk, quality 0.7 | 5.0% | 4.3% |
| Risk, quality 1.0 | 4.7% | 4.7% |

**Sturdy dogs helped more** (risk cut by 70 percent in the sturdier half and 30 percent in the frailer half; no dog's rate of aging changes). Share of trials that report slowed aging. Should be about 2.5 percent.

| Score | Without spread | With spread | Median estimated rate, with spread |
| --- | --- | --- | --- |
| No score | 52.2% | **87.2%** | 0.73 |
| Rank, quality 0.5 | 49.7% | 75.5% | 0.78 |
| Rank, quality 0.7 | 39.5% | **59.3%** | 0.82 |
| Rank, quality 0.9 | 21.7% | 27.5% | 0.90 |
| Rank, quality 1.0 | 9.2% | 10.0% | 0.95 |
| Risk, quality 0.7 | 45.3% | 64.8% | 0.81 |
| Risk, quality 1.0 | 12.3% | 14.2% | 0.94 |

**Risk cut by 50 percent in every dog** (no dog's rate of aging changes). Share of trials flagged. Should be about 5 percent.

| Score | Without spread | With spread |
| --- | --- | --- |
| No score | **24.2%**, all as faster aging | 5.0% |
| Rank, quality 0.7 | 11.5%, all as faster aging | 5.3% (check: 7.1%) |
| Risk, quality 0.7 | 10.5%, all as faster aging | 5.2% |

**Aging truly slowed 20 percent in every dog.** Share of trials that recognise it. Higher is better.

| Score | Without spread | With spread |
| --- | --- | --- |
| No score | 56.2% (check: 52.3%) | 69.7% |
| Rank, quality 0.7 | 62.2% | 67.5% |
| Rank, quality 1.0 | 72.8% | 73.3% |
| Risk, quality 0.7 | 63.5% | 68.8% |

## What this means

1. **The reworked formula is sound when nothing happens.** With no effect, the test with the spread setting flags 4 to 5 percent of trials at every score and quality. With no score it recovers the true spread: the median fitted value is 0.29 against the 0.3 used to make the data.
2. **The spread setting makes the unequal-benefit trap worse.** With sturdy dogs helped more, it raises false reports of slowed aging by about 35 points with no score, 26 points at quality 0.5 and 20 points at quality 0.7. The likely reason is the one given in experiment 017. The treatment itself makes treated dogs more unequal in risk than untreated dogs. The test allows one amount of spread for both groups. The extra spread among treated dogs has nowhere to go but into the estimated rate of aging.
3. **A very good score makes the two tests the same.** At quality 0.9 the gap is about 6 points, and at quality 1.0 there is none (about 10 percent either way). The score then leaves little spread over: the fitted setting sits at its lowest allowed value in a third or more of trials at quality 0.9 and more than half at quality 1.0.
4. **The spread setting is still needed for a plain cut in risk.** When every dog's risk is cut by 50 percent, the test without the spread setting reports faster aging in about 24 percent of trials with no score. A score of quality 0.7 halves that, to about 11 percent, and does not remove it. With the spread setting the rate is about 5 percent with no score and about 5 to 7 percent with a score of quality 0.7.
5. **The two mistakes pull in opposite directions.** Leave the spread setting out, and a plain risk cut is read as faster aging. Put it in, and an unequal benefit is read as slowed aging more often. A trial does not know beforehand which kind of treatment it has. So neither version can be called the safe default.
6. **The spread setting adds power for a true slowing, but that is the same property as point 2.** With aging truly slowed 20 percent, the test with the spread setting recognises it in about 68 to 70 percent of trials, against about 52 to 56 percent (no score) or 62 percent (quality 0.7) without it. The test with the spread setting is more sensitive to anything that looks like a change in slope, real or not. The gain is only a virtue if the benefit is known to be shared equally.
7. **The warning from experiment 034 still matters most.** A test of "does the benefit differ by score?" noticed the unequal benefit in 88 percent of trials at quality 0.7 (experiment 034). When that test says yes, a finding of slowed aging should not be relied on, with or without the spread setting.

For a real trial this suggests: adding a setting for left-over hidden differences to a score-adjusted slope test is not a repair for the unequal-benefit trap. Running both versions and comparing them may be more useful than choosing one, because they err in opposite directions. Whether that comparison can be turned into a sound test has not been tried here. This is a statement about the mathematics, not a comment on any trial's own analysis plan, which answers its own questions.

## Independent check

A separate agent wrote its own simulation and fitting code from a description of the model, with its own random draws and its own limits on the fitted settings. It ran about 1,500 trials for each of 21 figures before opening our code. All 21 were within 2 standard errors of ours. Its figures are in the tables above where they differ by 2 points or more, and in the predictions table.

It then read our code. It found no errors in the likelihood, its gradient, the simulation, the scores or the counting. It confirmed that every figure quoted matches `results.json`, and confirmed the two inputs (level 0.004967 a year and slope 0.3195 a year) against `dog_lifespan/results.json`. In its own runs the fit with a free rate was never worse than the fit without one.

It also checked the one change to experiment 012's formula. Dropping the term for frail dogs dying before entry is sound here: among dogs alive at any entry age, the hidden factor has the same relative spread, and entry age's own setting takes up the change in its average. It noted that this change in average is slightly curved with entry age and the setting is a straight line; no result suggests that matters.

## What the check changed

- **"Costs little power" became "adds power, which is the same property as the trap".** A first reading treated the gain in power (point 6) as a point in the spread setting's favour. The check pointed out that the gain and the extra false reports are one property seen twice. Point 6 now says so.
- **The spread test with a 50 percent cut and a score of quality 0.7 may be slightly high.** We found 5.3 percent; the check found 7.1 percent over 1,500 trials, and about 6.6 percent with both sets pooled, tilted towards faster aging. The table shows both. We describe it as "about 5 to 7 percent", not as exactly sound.
- **The gain in power with no score is larger in the check.** It found 52.3 percent without the spread setting against our 56.2. Both are shown.
- The check suggested why the spread setting adds power: hidden differences flatten the rise in risk in both groups by different amounts, and a test that ignores them blurs the contrast. It did not test this, and neither did we, so it is offered as a likely reason only.
- The check noted that with a score the fitted spread (about 0.1 at quality 0.7) is a working figure and not an estimate of anything in the model that made the data. We do not quote it as a finding.
- No claim was withdrawn.

## Comparison with published work

`papers/literature-check.md` names no published work on adjusting a slope test for a frailty score or for left-over hidden differences, and none was found in the searches run for it. It credits the idea of hidden frailty to Beard (1959) and Vaupel, Manton and Stallard (1979). No new search was run, because web pages could not be opened in this run. Raz et al. 2026 and McMillan et al. 2024 were tried again on 20 November 2026 and could not be opened, so experiments 022 and 023 still wait.

## Limits

1. Everything is assumed: the spread of hidden risk (0.3; nobody has measured it for dogs), the benefit sizes (70 and 30 percent, 50 percent, 20 percent slowing), the split into two halves, and the form of both scores.
2. Only the 10-year, 1,300-dog trial with one size of dog and no dropout was tested. Experiment 038 found the trap grows with the number of dogs; that was not repeated here.
3. Only sturdy dogs helped more was tested. Frail dogs helped more (experiment 037) was not.
4. The benefit may differ by score in a straight line only. The five score groups of experiment 034 were not combined with the spread setting.
5. The spread setting is one figure for both groups. A separate spread for treated dogs was not tried.
6. 600 trials per case: each figure is uncertain by 1 to 2 points. The same simulated dogs are used for every score within a case, so rows of a table are not independent.
7. When the score is very good the spread setting sits at its lowest allowed value in most fits. The check saw no sign that this upsets the 5 percent rate.

## Files

- `PREDICTION.md`: written and pushed before any code.
- `score_and_spread.py`: the experiment (about 3 minutes on 2 cores).
- `results.json`: every number above except the check's.
