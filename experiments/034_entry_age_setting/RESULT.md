# Experiment 034: result

**Giving entry age its own setting removes the fault found in experiment 032. With a treatment that does nothing and a perfect frailty score, the score-adjusted slope test flags about 5 percent of simulated trials, as it should, where the old version flagged about 11 percent. But the fault had been hiding part of the real problem. With the fault removed, a treatment that helps sturdy dogs more is read as slowed aging more often than experiment 032 reported: in about 40 percent of trials with a score of quality 0.7 (32 percent before), about 21 percent at quality 0.9 (8 to 12 before) and about 9 percent with a perfect score (2 to 3 before). Replacing the straight-line "benefit differs by score" term with five score groups helps only with a perfect score (about 5 percent). The entry-age setting also costs power when it is not needed.** All four predictions were confirmed.

This is a calculation under stated assumptions. It says how analyses behave on simulated dogs. It does not say that any treatment works, that any treatment helps sturdy dogs more, or that any real frailty measurement reaches any of these qualities. The spread of hidden risk (0.3), the form of the scores and every size of benefit are assumptions.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | No effect, rank score: the entry-age test flags 3 to 8% of trials at quality 0.9 and at 1.0, and its median estimated rate is between 0.98 and 1.03 | 4.8% and 5.5% (check: 6.0% and 4.7%); median rates 1.001 and 0.998 | Confirmed |
| P2 | Sturdy helped more, rank score: the entry-age test reports slowed aging in 8 to 22% at quality 1.0 and 28 to 45% at quality 0.7 | 9.0% (check: 9.3%); 39.7% (check: 41.5%) | Confirmed |
| P3 | Sturdy helped more, rank score: the five-group test reports slowed aging in 3 to 10% at quality 1.0, and is within 6 points of the straight-line version at quality 0.7 | 4.7% (check: 5.4%); 39.7% against 39.7% (check: 42.3% against 41.5%) | Confirmed |
| P4 | Aging truly slowed 20%, rank score: the entry-age test recognises it in 60 to 76% at quality 0.7 and 70 to 90% at quality 1.0; five groups no more than 5 points lower | 63.7% and 75.0% (check: 64.5% and 74.9%); five groups 61.5% and 72.8% | Confirmed |

## The three analyses

All three are the plain slope test with a frailty score included. Each asks: did the rate at which risk rises change in treated dogs?

- **Old**: experiment 032's version. One rate serves for both entry age and time in the trial. The benefit may differ by score in a straight line.
- **Entry-age**: the effect of entry age on risk is its own setting, separate from the rate at which risk rises during the trial.
- **Five groups**: as entry-age, with a separate risk cut for each fifth of the score in place of the straight line.

There are two kinds of score. The **rank score** is experiment 032's: the dog's rank in hidden risk among dogs of the same entry age, read with error. The **risk score** is the dog's hidden risk itself, read with error. Quality runs from 0 (useless) to 1 (perfect).

## Main tables

Each figure is the share of 600 simulated trials of 1,300 dogs followed for 10 years. Each is uncertain by about 1 to 2 percentage points.

**Treatment does nothing.** Share of trials flagged as a change in the rate of aging. Should be about 5 percent.

| Score | Old | Entry-age | Five groups |
| --- | --- | --- | --- |
| Rank, quality 0.5 | 5.7% | 5.0% | 4.5% |
| Rank, quality 0.7 | 4.8% | 5.5% | 5.0% |
| Rank, quality 0.9 | 7.8% (check: 9.9%) | 4.8% (check: 6.0%) | 5.2% |
| Rank, quality 1.0 | **10.7%** (check: 11.7%), nearly all as faster aging | 5.5% (check: 4.7%) | 5.2% |
| Risk, quality 0.7 | 5.5% | 4.8% | 4.8% |
| Risk, quality 1.0 | 4.7% | 4.7% | 5.0% |

**Sturdy dogs helped more** (risk cut by 70 percent in the sturdier half and 30 percent in the frailer half; no dog's rate of aging changes). Share of trials that report slowed aging. Should be about 2.5 percent.

| Score | Old | Entry-age | Five groups | Trial notices that the benefit differs by score (straight-line test) |
| --- | --- | --- | --- | --- |
| Rank, quality 0.5 | 46.2% | 49.7% | 47.8% | 53% |
| Rank, quality 0.7 | 29.8% | 39.7% (check: 41.5%) | 39.7% (check: 42.3%) | 88% |
| Rank, quality 0.9 | 9.7% | 21.5% (check: 21.1%) | 17.8% (check: 18.1%) | 99% |
| Rank, quality 1.0 | 1.8% (check: 2.6%) | 9.0% (check: 9.3%) | 4.7% (check: 5.4%) | 100% |
| Risk, quality 0.7 | 46.7% | 45.2% | 46.0% | 86% |
| Risk, quality 1.0 | 15.5% | 12.5% | 7.7% | 100% |

**Aging truly slowed 20 percent in every dog.** Share of trials that recognise it. Higher is better.

| Score | Old | Entry-age | Five groups |
| --- | --- | --- | --- |
| Rank, quality 0.7 | 68.8% | 63.7% | 61.5% |
| Rank, quality 1.0 | 64.3% | 75.0% | 72.8% |
| Risk, quality 0.7 | 79.0% | 64.0% (check: 67.0%) | 62.7% |
| Risk, quality 1.0 | 86.5% | 75.8% | 75.3% |

## What this means

1. **The explanation given in experiment 032 for the lean was right.** With a rank score, risk rises with entry age more slowly than with time in the trial once the score is known. Give entry age its own setting and the test with no effect goes back to about 5 percent at every score quality, with the estimated rate at 1.00.
2. **The lean had been flattering the repair.** The old test pushed every estimate towards "faster aging". When sturdy dogs are helped more, the error runs towards "slowed aging". The two partly cancelled. With the lean removed, false reports of slowed aging go up at every quality of the rank score: about 40 percent at 0.7, about 21 percent at 0.9, about 9 percent at 1.0. Experiment 032's statement that the repair works "at quality 0.9 or better" was too kind. With an honest analysis, even quality 0.9 leaves about one trial in five reporting slowed aging that is not there.
3. **Five score groups help only with a perfect score.** At quality 1.0 they bring false reports to about 5 percent (rank score) or 8 percent (risk score). At quality 0.9 they give about 18 percent, and at 0.7 they make no difference. The trouble at moderate quality is not the shape of the term. It is that the score does not know which dogs are frail.
4. **The entry-age setting has a price.** When the score is the risk score, the old test had no lean and did not need the extra setting. Adding it there takes power for a true 20 percent slowing from about 79 to about 64 percent at quality 0.7, and from about 87 to about 76 percent at quality 1.0. The reason: the test then learns the untreated rate only from what happens during the trial, not from how risk differs between dogs that entered at different ages. With the rank score at quality 0.7 the price is about 5 points; at quality 1.0 the entry-age test is the better one, because the old test's lean was working against it.
5. **A trial cannot know in advance which kind of score it has.** A real frailty measurement may behave like the rank score, like the risk score, or like neither. The entry-age setting is the safe choice against false alarms, and it costs 5 to 15 points of power.
6. **The warning still works, and is still the main use of a score.** The straight-line test of "does the benefit differ by score?" is answered yes in 88 percent of trials at quality 0.7, and in about 5 percent when the benefit is the same for every dog.
7. **Not predicted: the five-group version of that warning is faulty with a rank score.** With no treatment effect and a perfect rank score, a test of "are the five group cuts equal?" says no in 15 percent of trials, where 5 is expected. The check found the cause. Risk does not rise in a straight line with the rank score. The five group cuts let the treated dogs follow that bend while the untreated dogs are held to a straight line, and the test reads the bend as a difference in benefit. When the check gave both groups the same five score groups, the rate fell to about 6 percent. So a grouped test of this kind should group the score in both arms.

For a real trial this suggests: an analysis that adjusts a slope test for a frailty score should give entry age its own setting, and should expect to pay some power for it. Even then the score is a warning, not a repair. A finding of slowed aging in a trial where the benefit differs by the score should not be relied on, whatever the quality of the score. This is a statement about the mathematics, not a comment on any trial's own analysis plan, which answers its own questions.

## Independent check

A separate agent wrote its own simulation and fitting code from a description of the model, with its own random draws, wider limits on the fitted settings and three starting points for each fit. It ran ten cases with 1,500 trials each. Its figures are in the tables above. Every one of 28 comparisons was within 1.6 standard errors of ours. The largest gaps were the old test with no effect at quality 0.9 (9.9 percent against our 7.8) and the entry-age test for a true slowing with the risk score at quality 0.7 (67.0 against 64.0).

It then read our code and found no errors in the likelihood, the simulation, the scores, the five groups, the counting or the fitting limits. It confirmed that every figure quoted matches `results.json`, and confirmed the two inputs (level 0.004967 a year and slope 0.3195 a year) against `dog_lifespan/results.json`.

One departure from the plan, reported by the check itself: it opened our script after its own code was finished but while eight of its ten cases were still running. Its code was not changed afterwards.

## What the check changed

- Point 7 is the check's finding. We had the 15 percent figure; the check explained it and tested the explanation (about 6 percent with the score grouped in both arms; that figure is the check's own and is not a saved result).
- The headline was reworded. A first draft said the entry-age setting "repairs" the score-adjusted test. The check pointed out that it repairs only the false alarms with no effect; with sturdy dogs helped more it raises the false reports of slowed aging. The result now says so first.
- The check stressed that 5 percent of trials reporting slowed aging (five groups, perfect score) is still twice the 2.5 percent expected. That is stated in the table heading.
- Figures where the check differed by 2 points or more are shown with the check's value beside them.
- No claim was withdrawn.

## Comparison with published work

`papers/literature-check.md` names no published work on adjusting a slope test for a frailty score, and none was found in the searches run for it. It credits the idea of hidden frailty to Beard (1959) and Vaupel, Manton and Stallard (1979). That entry age and time in a study can carry different effects is a standard point in survival analysis; no source for it was opened in this run, so none is cited. No new search was run, because web pages could not be opened in this run (Raz et al. 2026 and McMillan et al. 2024 were tried again on 19 November 2026 and could not be opened).

## Limits

1. Everything is assumed: the spread of hidden risk (0.3; nobody has measured it for dogs), the benefit sizes (70 and 30 percent), the split into two halves, and the form of both scores.
2. Only the plain slope test was adjusted. The hidden-differences test of experiment 012 was left out, because its formula needs reworking when entry age has its own setting.
3. The true split between sturdy and frail dogs falls at the middle of the hidden-risk ranking, inside the middle one of the five score groups. A different true pattern, or a different number of groups, could give different figures.
4. Only the 10-year, 1,300-dog trial with one size and no dropout was tested.
5. The score is measured once, at entry, with plain bell-curve error.
6. 600 trials per case: each figure is uncertain by 1 to 2 points. The same simulated dogs are used for every score quality within a case, so rows of a table are not independent.
7. The effect of entry age is a straight line in the analysis. The true effect given the score is slightly curved; no harm from this was seen.

## Files

- `PREDICTION.md`: written and pushed before any code.
- `entry_age.py`: the experiment (about 2 to 5 minutes on 2 cores).
- `results.json`: every number above except the check's.
