# Experiment 016: result

**With their own published designs, the two running dog trials are well placed to see a large, lasting cut in death risk, and poorly placed to see a slower rate of aging. In the rapamycin trial as designed (580 dogs, one year of dosing, three years of follow-up), about 30 percent of placebo dogs die. A 20 percent slowing of aging during the dosing year shows up in about 6 percent of simulated trials, which is no better than chance. In a trial like Loyal's STAY study (1,300 dogs aged 10 or older, four years, treated throughout), about 83 percent of placebo dogs die. A lasting 20 percent cut in risk is detected in about 9 trials in 10, but a true 20 percent slowing of aging is recognised as a slope change in only about 3 in 10.** Two predictions were confirmed and two were half confirmed, half refuted.

This is a calculation under stated assumptions. It says what a design can detect on simulated dogs. It does not say whether any treatment works, and it does not question the trials' own power statements, which answer a different question (a large, lasting cut in risk).

**This run could not open the rapamycin trial's design paper** (Coleman et al. 2025). Its inputs are taken from the earlier reading recorded in this repository (`papers/literature-check.md`, `interventions/triad_published_inputs.py`). The independent check could not open it either.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | Rapamycin trial as published, aging slowed 20% in the dosing year: any-difference test under 10%, slope test under 10% | 6.4% and 4.7% (check: 6.0% and 5.1%) | Confirmed |
| P2 | Rapamycin trial as published, 50% risk cut in the dosing year that then fades: detected in 25 to 50% of trials; misread by the slope test in more than 20% | Detected in 16% (check: 16%). Misread in 38% (check: 37%) | **Half and half.** First part refuted: the trial sees a fading benefit less often than predicted. Second part confirmed |
| P3 | STAY-like trial: 60 to 85% of placebo dogs die; a 20% risk cut is detected in more than 85%; a 20% slowing is recognised in under 35% | 83% die under both guesses about entry ages. Risk cut: 92% and 89%. Slowing: 30% and 34.5% (check: 30% and 34%) | Confirmed. The last figure sits on the 35% line and is within noise of it |
| P4 | Rapamycin trial's dogs, treated throughout for 10 years: slope test recognises a 20% slowing in 80 to 95% of 1,300-dog trials and under 65% of 580-dog trials | 1,300 dogs: 78% (check: 79%). 580 dogs: 46% (check: 46%) | **Half and half.** First part refuted, narrowly: 1,300 dogs fall just short of 80%. Second part confirmed |

## Part A: the rapamycin trial as published

580 dogs, half treated for 12 months, all followed for 36 months. 1,000 simulated trials per row. About 30 percent of placebo dogs die in every row. Each figure is uncertain by about 1 to 1.5 percentage points.

| True effect | Any difference detected | Slope test reports a change | Two-period test reports a lasting benefit |
| --- | --- | --- | --- |
| None | 5.4% | 4.9% | 5.0% |
| 20% risk cut in the dosing year, then fades | 5.9% | 9.2% | 4.4% |
| 50% risk cut in the dosing year, then fades | 16.2% | **37.8%** | 4.7% |
| 20% risk cut that lasts all 3 years | 28.5% | 4.3% | 24.3% |
| 38% risk cut that lasts all 3 years | 78.7% | 5.2% | 71.5% |
| Aging slowed 20% in the dosing year, gain kept | 6.4% | 4.7% | 7.3% |
| Aging stopped in the dosing year, gain kept | 33.4% | 8.2% | 36.1% |

What this shows:

1. **The trial as designed can see a large lasting cut in risk.** A 38 percent cut is detected about 79 percent of the time, in line with the reconciliation in `interventions/triad_published_inputs.py` (81 percent).
2. **A benefit that ends when dosing ends is hard to see at all.** Even a halving of risk for the dosing year is detected in only 16 percent of trials. Few dogs die in the first year, so there is little to count. We predicted 25 to 50 percent; that was too hopeful.
3. **One year of modestly slowed aging is invisible** (6 percent, against 5 percent by chance). Even one year of aging stopped completely is detected in only about one trial in three.
4. **The slope test must not be used on this design.** It assumes the effect is constant from enrolment. When a 50 percent benefit fades, it reports a slope change in 38 percent of trials, and every one of those points to *faster* aging in treated dogs. These are errors of a wrong model, not detections. The two-period test of experiment 006 stays near 5 percent in the same cases.
5. **Early looks at the data cost almost nothing.** With looks after 200 and 400 dogs and a cautious stopping rule (ASSUMPTION, see limits), the 38 percent cut is detected in 77.5 percent of trials against 78.7 percent without looks, and 46 percent of such trials would stop early. Early looks do not add power; they only allow an earlier end.

## Part B: a trial like STAY

1,300 dogs aged 10 or older, half treated throughout, 4 years. 1,000 simulated trials per row. About 83 percent of placebo dogs die under both guesses about entry ages.

| True effect | Entry ages evenly 10 to 14: any difference / slope change | Entry ages shaped like the population aged 10 to 16: any difference / slope change |
| --- | --- | --- |
| None | 4.2% / 4.7% | 6.1% / 5.8% |
| Risk cut 20% | 92.2% / 5.6% | 88.6% / 4.0% |
| Risk cut 30% | 99.9% / 5.8% | 99.9% / 4.9% |
| Aging slowed 20% | 40.0% / 30.1% | 37.0% / 34.5% |
| Aging slowed 30% | 72.5% / 56.2% | 68.3% / 65.0% |

What this shows:

1. **A trial of this shape is strong for a lasting risk cut.** Most dogs die within the trial, so even a 20 percent cut is detected about 9 times in 10.
2. **It is weak for telling slowed aging apart.** A true 20 percent slowing is noticed as "some difference" in about 4 trials in 10, and recognised as a change in slope in about 3 in 10. This matches experiment 001 (28 percent). The guess about entry ages moves the answer by only a few points.
3. **Side calculation: what one more year of life takes.** A news report says STAY is designed to detect at least a one-year difference in lifespan. Under the UK death rates, a dog entering at age 10 to 14 gains one year on average from a lasting 46 percent cut in risk (risk ratio 0.54). A 20 percent cut adds about 0.34 years, about four months. This is our arithmetic, not a statement about STAY's design, which is not public.

## Part C: experiment 010 with the rapamycin trial's dogs

The same dogs as Part A, but treated throughout and followed for 10 years, when 98 percent of placebo dogs have died. No trial plans this; it is a what-if.

| Dogs | Slope test recognises a 20% slowing | False alarms (no effect / 50% risk cut) |
| --- | --- | --- |
| 580 | 45.9% | 4.9% / 5.4% |
| 1,300 | 77.7% (check: 79.4%) | 5.1% / 5.3% |

With this project's own inputs (UK death rates, entry ages 7 to 12), 1,300 dogs reached about 88 percent (experiments 010 and 012). With the rapamycin trial's published inputs (slope 0.29, entry ages 7 to 10) they reach about 78 percent. **So the trial size from experiment 010 depends on which death-rate inputs are used.** "About 1,300 dogs" is the low end; with these inputs somewhat more than 1,300 are needed for 80 percent. We did not test which input causes the drop (the lower slope or the narrower entry ages).

## Independent check

A separate agent wrote its own simulation and tests from a description of the model, before reading our code. It ran 2,000 to 3,000 trials per case with its own seeds and reproduced more than 30 numbers. All agreed within simulation noise; the largest gaps were about 1.6 standard errors. It checked the share of placebo dogs dying by direct integration (29.7, 82.9, 82.8 and 98.0 percent). It ran our simulator on 2 million dogs against exact formulas, and ran our fitting code on 1,600 of its own simulated trials: every decision matched its own. It found no errors and no departure from the method fixed in `PREDICTION.md`.

What the check changed or added:

- **Fitted values sometimes sit on a limit of the search range** (the slope ratio in 8 percent of the fading-benefit trials). The check refitted with wider limits and got the same decisions, so the 38 percent is not an artefact.
- **The published inputs give about 30 percent of placebo dogs dying in 3 years**, where the design paper is recorded as saying "approximately one-third". We now state this plainly. It is a rounding-size gap, not an error by anyone.
- **The second early look "after 400 dogs" is recorded only once in this repository** (in the literature check's list of changes). The first look, after 200 dogs, is recorded in its text. Neither could be re-read at the source.
- **"Treated throughout" for STAY is how the literature check reads the sources ("appears to dose continuously")**, not a confirmed fact. It is now listed as such.
- It added limits 1, 3, 5 and 7 below.

**Not checked at the source:** every fact about the two trials. The check could only confirm that our inputs match what this repository records from earlier readings.

## How this compares with published work

`papers/literature-check.md` found no published power analysis for detecting a slope change in a companion-dog trial in the searches run, so there is no published figure to compare the slope results with. The rapamycin trial's design paper (Coleman et al. 2025) states greater than 80 percent power for a lasting cut in the level of risk; our 79 percent for a lasting 38 percent cut is consistent with the reading in `papers/literature-check.md`, and that reading is an interpretation, as recorded there. That a benefit can vanish, partly last, or rebound after treatment stops is published in animal studies (Mair et al. 2003; Hahn et al. 2019; McCracken et al. 2020, as cited in the literature check; not opened in this run). This experiment simulates only vanishing and fully lasting benefits.

## Limits

1. **Best case throughout.** The simulated dogs follow the doubling law exactly, nobody drops out, every owner gives the treatment, and the analysis has the right form except where stated. Real detection rates would be lower.
2. **ASSUMPTION: the any-difference test** is a two-sided log-rank test at the 5 percent level. The design paper is recorded as naming a Cox model; the two are close for a single comparison.
3. **ASSUMPTION: the early-look rule.** The stopping rule is ours (O'Brien-Fleming type), and each look uses complete three-year records of the first 200 and 400 dogs. Real looks happen by calendar date with partial records.
4. **ASSUMPTIONS in Part B:** death rates from the UK table of all dogs (trial dogs are screened and may die more slowly); no dropout; both sets of entry ages; treatment throughout. Published, as recorded in the literature check: 1,300 dogs, age 10 or older, about four years.
5. **The UK doubling law overstates death rates from about age 16** (experiment 004). Many STAY-like dogs reach those ages, so the 83 percent is probably somewhat high.
6. **The rapamycin trial has enrolled far fewer than 580 dogs** at the last published count (210 in September 2025). Part A describes the design, not the trial as it stands.
7. **The two-period test can be biased by which dogs survive the dosing year** (Bartlett et al. 2020, as cited in the literature check). No sign of it appears here, but a rebound after stopping was not simulated.
8. Part C is a what-if with one effect size and two trial sizes.

## Inputs

- Rapamycin trial: level 0.0066, slope 0.29 per year, entry ages 7 to 8 (half), 8 to 9 (quarter), 9 to 10 (quarter), 580 dogs, 12 months of dosing in 36: Coleman et al. 2025, GeroScience 47:2851, as recorded in `interventions/triad_published_inputs.py` and `papers/literature-check.md`. Not re-opened.
- STAY: 1,300 dogs, aged 10 or older, about 4 years: DVM360 (13 January 2026) and Loyal (9 October 2024), as recorded in `papers/literature-check.md`. Not re-opened.
- UK level and slope: `dog_lifespan/results.json`.
- Slope test: `experiments/001_trial_fingerprint/fingerprint.py`. Two-period test: `experiments/006_stop_and_watch/two_period.py`. Log-rank test and entry-age mix: `interventions/triad_published_inputs.py`.

## Files

`real_designs.py` (about 2 minutes), `results.json`, `PREDICTION.md`.
