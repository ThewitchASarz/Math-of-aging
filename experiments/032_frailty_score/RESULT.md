# Experiment 032: result

**A frailty score measured at entry helps only if it is very good. When a treatment helps sturdy dogs more than frail ones, the plain slope test reports slowed aging in 56 percent of simulated 10-year trials, although no dog aged more slowly. Adding the score and a "benefit differs by score" term brings that to 52 percent with a weak score (quality 0.3), 46 percent at quality 0.5, 32 percent at 0.7, about 8 to 12 percent at 0.9 and about 2 to 3 percent with a perfect score. But the same analysis has a fault of its own: with a very good score it leans towards "faster aging" when the treatment does nothing (10 percent of trials flagged where 5 is expected). The warning sign is easier to get than the repair: the trial notices that the benefit differs by score in 88 percent of trials at quality 0.7.** Two predictions were confirmed and two were half confirmed and half refuted.

This is a calculation under stated assumptions. It says how one analysis behaves on simulated dogs. It does not say that any treatment works, that any treatment helps sturdy dogs more, or that any real frailty measurement reaches any of these qualities. The spread of hidden risk (0.3), the form of the score and every size of benefit are assumptions.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | Sturdy helped more, no score: plain test flags slowed aging in 50 to 65% of trials, hidden-differences test in 82 to 93% | Plain 56.2% (check: about 58%); hidden-differences 87.8% (check: 87 to 91%) | Confirmed |
| P2 | Score terms in the plain test: still more than 30% at quality 0.5, and 20 to 45% at quality 0.7 | 46.5% (check: 44.7%); 32.2% (check: 35.3%) | Confirmed |
| P3 | Score terms in the plain test: 8 to 25% at quality 0.9, and 3 to 10% at quality 1.0 | 8.0% (check: 12.0%); 1.8% (check: 3.0%) | **Half and half.** Quality 0.9 is as predicted. A perfect score leaves fewer false reports than predicted, but part of that is the analysis leaning the other way (see point 3 below) |
| P4 | At quality 0.7: no effect, both tests with score terms flag 3 to 8%; aging truly slowed 20%, each test with score terms is no more than 5 points below the same test without | No effect: 4.7% and 5.5%. Slowing: plain 70.5% against 72.0% (1.5 points lower); hidden-differences 75.8% against 82.3% (6.5 points lower) | **Half and half.** No extra false alarms. The plain test loses little in our run; the hidden-differences test loses about 6 points, more than predicted. The check found losses of about 7 and about 10 points |

## Main table

Share of 600 simulated trials in which each test reports slowed aging. Every trial has 1,300 dogs, half treated from enrolment, followed for 10 years, with widely spread hidden risk (0.3, ASSUMPTION). Treated dogs in the sturdier half have their risk cut by 70 percent and the frailer half by 30 percent (ASSUMPTION). No dog's rate of aging changes, so every figure should be near 2.5 percent. Each figure is uncertain by about 1 to 2 percentage points.

| Quality of the score (0 = useless, 1 = perfect) | Plain test with score terms | Hidden-differences test with score terms | Trial notices that the benefit differs by score |
| --- | --- | --- | --- |
| no score used | 56.2% | 87.8% | cannot |
| 0.3 | 52.3% | 85.3% | 20% |
| 0.5 | 46.5% | 74.2% | 50% |
| 0.7 | 32.2% | 49.5% | 88% |
| 0.9 | 8.0% (check: 12.0%) | 14.0% (check: 17%) | 100% |
| 1.0 | 1.8% | 2.7% | 100% |

Other cases, with a score of quality 0.7 unless stated. Share of trials flagged as a change in the rate of aging, in either direction.

| True effect | Plain | Hidden-differences | Plain with score terms | Hidden-differences with score terms | Should be |
| --- | --- | --- | --- | --- | --- |
| None | 5.2% | 5.3% | 4.7% | 5.5% | about 5% |
| None, score quality 0.9 | | | 6.5% | 6.7% | about 5% |
| None, score quality 1.0 | | | **10.0%**, nearly all as faster aging | **10.5%** | about 5% |
| 50% risk cut in every dog | 23.2% | 5.7% | 17.5% | 7.5% | about 5% |
| Frail dogs helped more (cuts of 70% and 30%) | 77.2%, all as faster aging | 30.3% | 32.0% (check: 39%) | 15.5% (check: 20%) | about 5% |
| Aging slowed 20% in every dog | 72.0% | 82.3% | 70.5% | 75.8% | high |
| Aging slowed 20%, score quality 1.0 | | | 64.2% | 65.7% | high |

The rows with score quality 0.9 and 1.0 for "none" and 1.0 for "aging slowed" were added after the first run, when the estimated rate with no effect was seen to sit a little above 1. They were not in the prediction.

## What this means

1. **A middling frailty score does not protect a slope test.** A score of quality 0.5 leaves nearly half of trials reporting slowed aging that is not there, and 0.7 leaves about a third. The false reports fall below about 10 percent only at quality 0.9 or better. Quality here means how closely the score tracks the dog's true standing in hidden risk among dogs of its entry age. Nobody has measured what quality a real frailty index, a vet's score or a blood test reaches for dogs; this experiment cannot say.
2. **The warning is cheaper than the repair.** The simple question "does the benefit differ by score?" is answered yes in 50 percent of trials at quality 0.5 and 88 percent at quality 0.7, and in about 5 to 7 percent when the benefit is in fact the same for every dog. So a moderately good score can tell a trial that a finding of slowed aging is not safe to trust, even though it cannot make the finding safe.
3. **The repair has a fault of its own.** With no treatment effect, the test with score terms flags 4.7 percent of trials at quality 0.7, 6.5 percent at 0.9 and 10.0 percent at 1.0, nearly all as faster aging. The reason is in how the score was built. It ranks a dog among dogs of the same entry age. A dog of a given rank that enters older has lower hidden risk, because frail dogs are rarer among older entrants. So, given the score, risk rises with entry age more slowly than it rises with time in the trial. The analysis gives untreated dogs one rate for both and treated dogs a free rate for time in the trial, and the free rate comes out too high. The check confirmed this in one trial of 52,000 dogs: estimated rate 1.066 times the untreated one with no effect, and 1.008 when entry age was given its own setting. This means the very low figure for a perfect score (1.8 percent) is partly two errors pulling against each other.
4. **With a score that reads hidden risk itself, the lean disappears but the repair is incomplete.** An extra run, not in the prediction, used the dog's hidden risk itself as the score, read without error (400 trials). With no effect the test flags 5.8 percent of trials and the estimated rate is 1.00. With sturdy dogs helped more it still reports slowed aging in 15 percent of trials (check: 14 percent). A straight-line "benefit differs by score" term does not fully describe a benefit that jumps between two halves. A true 20 percent slowing is recognised in 87 percent of these trials, more than without a score (72 to 82 percent), because the score explains differences between dogs.
5. **The score terms cost little when not needed, at moderate quality.** At quality 0.7 there are no extra false alarms and the plain test loses about 1.5 to 7 points of power for a true slowing. The hidden-differences test loses about 6 to 10 points. At quality 1.0 with the rank-based score, power falls to about 64 percent, because of the lean described in point 3.
6. **The hidden-differences test stays the worse of the two when sturdy dogs are helped more**, at every score quality below 1 (for example 49.5 against 32.2 percent at quality 0.7). That repeats experiment 017.

For a real trial this suggests: measuring frailty at entry is worth planning, mainly as a warning. If the benefit differs by the score, a slope-test finding of slowed aging should not be relied on. If the benefit does not differ by a moderately good score, that is weak reassurance only, since at quality 0.5 half of trials miss a large difference. And an analysis that adjusts for a score should let risk depend on entry age separately from time in the trial. This is a statement about the mathematics, not a comment on any trial's own analysis plan, which answers its own questions.

## Independent check

A separate agent wrote its own simulation and fitting code from a description of the model before opening ours, with its own random draws and a different way of computing gradients, and ran every case (300 trials each). Its figures are in the tables above. All comparisons were within about two standard errors. The largest gaps were the plain test with score terms at quality 0.9 (12.0 percent against our 8.0) and the frail-helped-more case with score terms (39 against 32 percent).

It then read our code and found no errors in the likelihood, the simulation, the counting or the fitting. It confirmed that every figure quoted matches `results.json`. It confirmed the two inputs (level 0.004967 a year and slope 0.3195 a year) against `dog_lifespan/results.json`, the UK fit at ages 2 to 19.

It tested our explanation of the lean towards faster aging (point 3) and agreed with it, by reasoning and with the 52,000-dog trial described there. That was one large trial, not a repeated test.

## What the check changed

- The figure at quality 0.9 is stated as about 8 to 12 percent, not 8.
- The power cost of the score terms at quality 0.7 is stated as a range (plain test 1.5 to 7 points, hidden-differences test 6 to 10), because the check's run found more than ours.
- The frail-helped-more figures with score terms are stated with the check's higher values beside them.
- The explanation in point 3 was a suggestion before the check and is now supported by a direct test.
- No claim was withdrawn.

## Comparison with published work

`papers/literature-check.md` names no published work on using a measured frailty score to protect a slope test, and none was found in the searches run for it. It credits the idea of hidden frailty itself to Beard (1959) and Vaupel, Manton and Stallard (1979), as cited by Missov et al. 2015. This experiment uses that idea; it does not add to it. No new search was run for this experiment, because web pages could not be opened in this run.

## Limits

1. Everything is assumed: the spread of hidden risk (0.3; nobody has measured it for dogs), the benefit sizes (70 and 30 percent), the split into two halves, and the form of the score.
2. The score ranks dogs within their entry age. That choice causes the lean in point 3. A score built another way would behave differently, as point 4 shows.
3. The "benefit differs by score" term is a straight line. Other forms were not tested.
4. Only the 10-year, 1,300-dog trial with one size and no dropout was tested. Not a 4-year trial, not the designs of the trials now running.
5. The score is measured once, at entry, and its error is plain bell-curve noise unrelated to anything else.
6. 600 trials per case (400 in the extra run): each figure is uncertain by 1 to 2 points, and differences of a few points between analyses are not established.
7. The two rows at quality 0.9 and 1.0 with no effect, and the extra run, were added after seeing the first results. They are reported, not scored.

## Files

- `PREDICTION.md`: written and pushed before any code.
- `frailty_score.py`: the experiment (about 5 minutes on 2 cores).
- `results.json`: every number above.
