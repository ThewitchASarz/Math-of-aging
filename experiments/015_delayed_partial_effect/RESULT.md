# Experiment 015: result

**A slow start does not weaken the slope test when the treatment truly slows aging, but it can fool the test when the treatment only cuts risk. A 50 percent risk cut that takes 2 years to reach full strength is reported as slowed aging in about 23 percent of simulated 1,300-dog, 10-year trials. Leaving the first 2 years out of the analysis brings that back to about 5 to 7 percent, but then a true 20 percent slowing is recognised in only about 51 to 54 percent of trials, down from about 79 to 81. Owners who do not give the treatment cost power, though less than predicted: with 20 percent of treated dogs never getting it, a true slowing is recognised in about 67 to 69 percent of trials; with 10 percent of owners stopping each year, in about 46 to 47 percent.** Two predictions were confirmed in full, one was confirmed on the edge, and one was half confirmed and half refuted.

This is a design calculation under stated assumptions. It shows what the mathematics and this repository's fitted inputs imply. It does not say that any treatment works, and it does not replace a real trial. Every delay and every share of owners here is an assumption, not a measurement.

## Predictions scored

"Check" is the independent check's own figure, from its own code.

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | A 20% slowing that builds up over 1 year is recognised in 72 to 82% of trials; over 2 years in 66 to 80% | 78% (check: 81%); 79% (check: 80%). With no delay: 79% (check: 81%) | Confirmed |
| P2a | 10% never given: 64 to 76%. 10% of owners stopping each year: 45 to 65% | 74% (check: 72%); 46% (check: 47%) | Confirmed |
| P2b | 20% never given: 52 to 66%. 30% never given: 40 to 56%. 1-year build-up and 20% never given, 1,300 dogs: 50 to 65% | 67% (check: 69%); 61% (check: 60%); 72% (check: 69%) | **Refuted**: owners who never treat cost less than predicted |
| P2c | Combined case: 1,600 dogs under 78%; the smallest trial on the grid reaching 80% has 2,000 or 2,600 dogs | 76% (check: 78%); 2,000 dogs reach 84% (check: 84%) | Confirmed |
| P3 | A 50% risk cut building up over 2 years is reported as a slope change in 20 to 50% of trials, middle estimate of the rate multiplier 0.90 to 0.97; over 1 year in 8 to 25%; with 20 or 30% never given (no delay) in under 12% | 23% (check: 22%), multiplier 0.91; 8.6% (check: 8.6%); 7% and 8% | Confirmed, the 1-year figure on the edge of its band |
| P4 | Analysis from year 2: the 2-year build-up is flagged in 3 to 8% of trials, no effect in 3 to 8%, and an immediate 20% slowing is recognised in 45 to 68% | 5.5% (check: 6.9%); 6.2% (check: 5.5%); 51% (check: 54%) | Confirmed |

P2 as a whole is half confirmed and half refuted.

## Main tables

All trials: 1,300 dogs unless stated, 10 years, entry ages 7 to 12, mix of sizes, hidden differences between dogs, 5 percent of dogs lost each year. 800 simulated trials per case here. Each figure is uncertain by about 1.5 percentage points, so cases that differ by less than about 4 points should not be ranked.

**A true 20 percent slowing of aging: how often the slope test recognises it**

| What happens | Recognised | Check | Middle estimate of the rate multiplier (truth: 0.80) |
| --- | --- | --- | --- |
| Works at once, every owner treats | 79% | 81% | 0.79 |
| Builds up over 1 year | 78% | 81% | 0.80 |
| Builds up over 2 years | 79% | 80% | 0.80 |
| 10% of treated dogs never get it | 74% | 72% | 0.81 |
| 20% never get it | 67% | 69% | 0.82 |
| 30% never get it | 61% | 60% | 0.84 |
| 10% of owners stop each year | 46% | 47% | 0.87 |
| Builds up over 1 year, 20% never get it | 72% | 69% | 0.82 |

**The combined case (1-year build-up, 20 percent never given): how many dogs**

| Dogs | Recognised | Check |
| --- | --- | --- |
| 1,300 | 72% | 69% |
| 1,600 | 76% | 78% |
| 2,000 | 84% | 84% |
| 2,600 | 92% | 93% |

**No true change in the rate of aging: how often the slope test reports one anyway**

| What happens | Usual analysis | Middle estimate of the rate multiplier (truth: 1.00) | Analysis from year 2 |
| --- | --- | --- | --- |
| No effect | 5% (check: 7%) | 1.00 | 6% (check: 6%) |
| 50% risk cut, at once | 6% (check: 5%) | 1.02 | not run |
| 50% risk cut, builds up over 1 year | 9% (check: 9%) | 0.96 | 7% (check: 7%) |
| 50% risk cut, builds up over 2 years | **23%** (check: 22%) | 0.91 | 5.5% (check: 7%) |
| 50% risk cut at once, 20% never given | 7% (check: 6%) | 0.98 | not run |
| 50% risk cut at once, 30% never given | 8% (check: 8%) | 0.96 | not run |

**The price of the repair.** With the analysis started at year 2, a true 20 percent slowing that works at once is recognised in 51 percent of trials (check: 54 percent), against 79 (check: 81) with the usual analysis. About 72 percent of the dogs are still being followed at year 2.

## What this means

1. **A delay in a true slowing is harmless to the test.** Most of the test's information comes from the later years, when the two groups' risks have drifted well apart. Whether the drift starts on day one or builds over a year or two makes no measurable difference. The trial sizes of experiments 010 and 011 do not need to grow for this.
2. **A delay in a plain risk cut is a trap.** While the cut is building up, the treated group's risk falls further below the untreated group's each month. The slope test reads that widening gap as a slower rise of risk. With a 2-year build-up it reports slowed aging in about one trial in four when the rate of aging has not changed at all. A 1-year build-up does much less harm (about 9 percent).
3. **Leaving out the build-up period repairs it, and costs about a third of the power.** A trial would need to be larger to afford this: the from-year-2 analysis with 1,300 dogs has roughly the power of a coin toss.
4. **Owners who never treat dilute the result but do not distort it much.** Each 10 percent of treated dogs that never get the treatment costs about 6 points of power, less than our simple dilution argument predicted. The estimate of the slowing shrinks toward "no effect" (0.82 where the truth is 0.80 among dogs that are treated, with 20 percent never given). With a plain risk cut, owners who never treat lift false alarms only slightly (7 to 8 percent).
5. **Owners who stop along the way cost the most.** Ten percent stopping each year takes the test from about 80 percent to about 46 percent. Over 10 years most owners have stopped, and the analysis still counts their dogs as treated.
6. **With a 1-year build-up and 20 percent never given, the trial needs between 1,600 and 2,000 dogs**, against about 1,300 in experiment 011. The grid is coarse; 2,000 is the smallest size tested that reaches 80 percent.

For a real trial that plans a slope test this means: decide in advance how long the treatment is expected to take to reach full strength, and plan a second analysis that leaves that period out; record whether and when owners stop giving the treatment; and size the trial for the owners who will stop.

## How this compares with published work

`papers/literature-check.md` names no published work on how a slow start or owners not treating affects a test of the rate of aging in a dog trial, and found no published figure for how often such a trial would recognise a slower rise of risk. That is a statement about a limited search. The general ideas used here, analysing dogs in the group they were assigned to and leaving out an early period, are standard practice in clinical trials and are not ours. Neither running trial has published a plan for a slope test, so nothing here describes what those trials will do; their own power statements answer a different question.

## Independent check

A separate agent wrote its own simulation and fitting code from a description of the model, before reading ours. It ran every case, with 800 to 3,000 trials each. All figures agreed within simulation noise. It then read our code and found no errors in the simulation, the solving of death times on a grid, the formula for owners who stop, the from-year-2 analysis, the likelihood or the bookkeeping. It matched all four inputs exactly against the files in this repository.

What the check changed:

- **A build-up is not said to help.** Our run showed the combined case (72 percent) above "20 percent never given" alone (67 percent), as if the delay helped. The check ran both with 3,000 trials and got 68.6 and 68.8 percent. The gap was noise. Likewise our 78 percent for a 1-year build-up against the check's 81 is one low draw.
- **Small differences are not ranked.** At 800 trials a figure is uncertain by about 1.5 points. The check's own first run of the stopping case gave 52 percent and its 3,000-trial re-run gave 47, close to our 46.
- **The trial size for the combined case is stated as a range**, 1,600 to 2,000 dogs.
- It ran the from-year-2 analysis for a slowing that builds up over 1 or 2 years (55 and 57 percent recognised). Those figures are the check's own and are not saved results.
- It added limits 2, 4 and 6 below.
- It noted that our script tests its grid method against the exact formula only for the case with no build-up. The check's separate grid gave the same answers for the build-up cases.

**Not checked against the outside source:** the size table's paper (Montoya et al. 2023) was not re-opened; the shares were matched against this repository's file.

## Limits

1. **Everything about delay and owners is an ASSUMPTION:** the straight-line build-up, its length (1 or 2 years), the shares never given (10, 20, 30 percent), the stopping rate (10 percent a year), and that a dog keeps the gain already made when its owner stops. No published figure was opened for any of them.
2. **The simulated effects have exactly the shape the analysis expects** once they are at full strength (a doubling law per size group; a pure slowing or a pure risk ratio). Effects of other shapes may behave differently, and a real analysis will fit less well.
3. **The misreading in point 2 is about the kind of effect, not its existence.** In those trials the treatment does cut deaths. What is wrong is the label "slowed aging".
4. **The from-year-2 repair was tested only here.** It keeps only dogs still followed at year 2, and those are the sturdier dogs in each group. In these simulations that did not disturb the test (5 to 7 percent false alarms), but it was tested for one size of hidden differences only. It also needs the build-up time to be known in advance; a cut-off shorter than the true build-up was not tested.
5. **Owners who stop were combined only with an immediate slowing.** Stopping with a risk cut, stopping linked to the dog's health, and owners of untreated dogs obtaining the treatment were not simulated.
6. **The combined case was run on a coarse grid of sizes** (1,300, 1,600, 2,000, 2,600).
7. All the assumptions of experiment 011 carry over: dropout of 5 percent a year unrelated to health, entry ages 7 to 12, size mix equal to the US population mix, the hidden risk factor of experiment 004 applied inside every size class, size-class rates that are a calibration.
8. **Sources not opened in this run:** Raz et al. 2026 and McMillan et al. 2024 were tried on 18 November 2026 and could not be opened (the permission request was not answered, as the run was unattended). Experiments 022 and 023 still wait.

## Inputs

- Size-class level and rates of aging: `dog_lifespan/data/size_class_gompertz_calibration.csv` (this repository's calibration).
- Size mix: `dog_lifespan/data/montoya2023_size_e0.csv` (Montoya et al. 2023, Frontiers in Veterinary Science 10:1082102).
- Spread of the hidden risk factor, 0.102: `experiments/004_old_dogs/results.json`.
- Trial scenario and slope test: `experiments/011_realistic_trial/realistic_trial.py`; likelihood and fitting: `experiments/007_breed_differences/breeds.py`.

## Files

`delayed_partial.py` (about 4.5 minutes), `results.json`, `PREDICTION.md`.
