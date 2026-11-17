# Experiment 019: result

**If death risk rebounded above the untreated level after dosing stopped, a trial that doses for one year and watches for three would tend to read a treatment that helped while it was given as harmful, or miss the pattern altogether. With risk halved in the dosing year and then 1.5 times the untreated level, the trial's main test fires in about 28 percent of simulated trials, and every time it fires the treated dogs are the ones that died more. The two-period test of experiment 006 sees that lasting rebound in about 76 percent of trials. It sees a short rebound (double risk for six months) in only about 30 percent, and the main test does not see it at all (5.5 percent, which is chance). When half of a benefit lasts, the main test fires in about 63 percent of trials but the benefit is recognised as lasting in only about 35 percent. The worry that the surviving dogs are a different mix turns out to be too small to matter in this design.** All four predictions were confirmed.

This is a calculation under stated assumptions. No rebound has been reported in dogs in the sources this project has read, and nothing here says any treatment rebounds or works. Every effect size is an assumption, chosen to see what the design could detect. The results do not question the trial's own power statement, which answers a different question (a large, lasting cut in risk).

**This run could not open the rapamycin trial's design paper** (Coleman et al. 2025). Its inputs are taken from the earlier reading recorded in this repository. The independent check could not open it either.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | Lasting rebound (risk halved in the dosing year, then 1.5 times untreated): the any-difference test fires in 20 to 45% of trials, and in at least 90% of those the treated group died more | Fires in 28.0% (check: 26.5%). Treated group died more in 100% of those (check: 99.8%) | Confirmed |
| P2 | The two-period test reports "worse after stopping" in 60 to 85% of trials for the lasting rebound, and in under 35% for the short rebound | 75.9% (check: 76.4%) and 30.3% (check: 29.5%) | Confirmed |
| P3 | Half the benefit lasts (risk halved, then cut by a quarter): any-difference test fires in 50 to 75%; two-period test reports a lasting benefit in 20 to 45% | 62.5% (check: 63.7%) and 35.2% (check: 37.1%) | Confirmed |
| P4 | Benefit fades, dogs differ in hidden ways (spread 0.30): two-period test reports an after-dosing difference in under 8% of trials; the calculated after-dosing risk ratio between survivors is below 1.03 | 5.8% (check, with 8,000 trials: 5.1%). Ratio 1.013 comparing dogs of the same entry age, 1.015 comparing the whole groups | Confirmed |

## The table

580 dogs, half treated for 12 months, all followed for 36 months. 1,000 simulated trials per row. About 30 percent of placebo dogs die in three years; about 8 percent of placebo dogs die in the dosing year. Each figure is uncertain by about 1 to 1.5 percentage points. "Risk multiplier" is the treated dogs' risk as a multiple of untreated risk at the same age.

| True effect (risk multiplier: dosing year / months 12 to 18 / months 18 to 36) | Any difference detected (share where treated died more) | Two-period test: worse after stopping | Two-period test: benefit lasts | Typical after-dosing ratio estimated |
| --- | --- | --- | --- | --- |
| No effect (1 / 1 / 1) | 4.6% (2.4%) | 2.7% | 2.1% | 1.00 |
| Fades (0.5 / 1 / 1) | 16.3% (0.4%) | 2.5% | 2.0% | 0.99 |
| Lasting rebound (0.5 / 1.5 / 1.5) | 28.0% (28.0%) | **75.9%** | 0.0% | 1.51 |
| Short rebound (0.5 / 2 / 1) | 5.5% (3.7%) | 30.3% | 0.1% | 1.26 |
| Small lasting rebound (0.8 / 1.25 / 1.25) | 11.1% (10.9%) | 26.7% | 0.0% | 1.25 |
| Rebound with no benefit first (1 / 1.5 / 1.5) | 55.6% (55.6%) | 75.3% | 0.0% | 1.50 |
| Half the benefit lasts (0.5 / 0.75 / 0.75) | 62.5% (0.0%) | 0.0% | 35.2% | 0.75 |
| Half of a smaller benefit lasts (0.7 / 0.85 / 0.85) | 26.0% (0.0%) | 0.0% | 15.3% | 0.85 |

With hidden differences between dogs (the treatment truly does nothing after the dosing year):

| True effect | Spread of hidden differences | Any difference detected | Two-period test reports an after-dosing difference (worse / lasting) |
| --- | --- | --- | --- |
| No effect | 0.10 | 5.9% | 4.9% (2.4% / 2.5%) |
| Fades (0.5 / 1 / 1) | 0.10 | 13.2% | 5.2% (3.1% / 2.1%) |
| No effect | 0.30 | 4.9% | 5.7% (2.8% / 2.9%) |
| Fades (0.5 / 1 / 1) | 0.30 | 14.3% | 5.8% (3.6% / 2.2%) |

## What this shows

1. **A rebound can turn a real benefit into an apparent harm.** Few dogs die in the dosing year (about 8 percent), and many more die in the two years after. So the two years after carry most of the weight in a test over all three years. A halving of risk while dosing, followed by a lasting 1.5-fold rebound, leaves more treated dogs dead by the end (36 percent against 30) and the main test, when it fires, always points to harm. Nothing in that one result would show that the treatment helped while it was given.
2. **A benefit while dosing can hide a rebound, and a rebound can hide a benefit.** With the same lasting rebound and no benefit first, the main test fires in 56 percent of trials. With the benefit first, 28 percent. With a short rebound, the two cancel almost exactly and the main test fires at the chance rate.
3. **The two-period test is the right tool, with limited reach.** It separates the dosing year from the time after, so it finds the lasting 1.5-fold rebound in about three trials in four, with or without a benefit first, and it almost never called a rebound a lasting benefit (1 trial in 1,000 in one row) or the reverse. But a smaller lasting rebound (1.25-fold) is found in about one trial in four, and a short one in about three in ten.
4. **One after-dosing number cannot tell a short, large rebound from a lasting, small one.** The test fits a single ratio for the whole period after dosing. Under the short rebound (double for six months, then normal) it typically reports about 1.26, which is neither true value. A lasting 1.25-fold rebound gives the same number. Telling them apart needs the after-dosing period split further, which this design has too few deaths for.
5. **A benefit that half lasts is usually not recognised as lasting.** The main test sees it in about 63 percent of trials, mostly because of the lasting part. The two-period test confirms that something lasts in about 35 percent. For a smaller benefit (30 percent cut, half lasting) the figures are 26 and 15 percent. So "no lasting benefit found" in a trial of this size is weak evidence that nothing lasted.
6. **Survivor mix is not a practical problem here.** If treatment keeps frail dogs alive through the dosing year, the treated survivors are slightly frailer. In this design the effect is tiny, because only about 8 percent of dogs die in the dosing year. Even with widely spread hidden differences (0.30), treated survivors die at 1.013 times the placebo rate afterwards when dogs of the same entry age are compared, and 1.015 times when the whole groups are compared (the second figure also counts a small shift in the ages of the survivors). The two-period test's false alarms stay at about 5 percent. In an older group (entry at 10 to 14, UK death rates, an assumption) the whole-group ratio would be about 1.05, so the issue grows where more dogs die during dosing.
7. **The slope test again misreads all of these.** It reports a slope change in 77 percent of lasting-rebound trials, always as faster aging. This repeats the warning of experiments 003 and 016 and was not scored.

## Independent check

A separate agent wrote its own simulation, log-rank test and two-period fit from a description of the model, before reading our code. It ran 8,000 trials per case with its own seeds. All 30 or so numbers it compared agreed within simulation noise (largest gap 1.4 standard errors). It tested our simulator on 2 million dogs against exact formulas, with and without hidden differences, reran our script and got identical output, and compared our fitting code with its own on 2,000 trials: no decision differed. It found no coding error and no departure from `PREDICTION.md`.

What the check changed or added:

- **The survivor ratio had two meanings mixed together.** Our calculation compared dogs of the same entry age (1.013). Our brute-force check on 2 million dogs measured the whole groups (1.015). These are different quantities, so one was not a check of the other. We added the whole-group calculation to the code (1.0153, against 1.0150 by brute force) and now report both. P4 holds under either.
- **The 5.8 percent in the "fades" world with hidden differences is noise around 5 percent.** With 8,000 trials the check got 5.1 percent. It found a very small tilt towards "worse after stopping" verdicts in the fading worlds (about 2.8 percent against 2.5 by chance), present even without hidden differences, so it is not survivor mix.
- **A small flaw in how we weighted entry ages in the side calculation** (0.4988 for the youngest group in place of 0.5) was fixed. No reported figure changed at the precision shown.
- **In the short-rebound world our fit sometimes stops slightly early** and about 3 percent of fits touch a search limit. No decision changed in the check's comparison.
- It supplied points 4 and the definition in point 6 above, and limit 3 below.

**Not checked at the source:** the trial's design inputs and the three papers on what happens after stopping. The check's request to open the design paper was not answered, and it did not work around that. It confirmed only that our inputs match what this repository records. The hidden-differences figure of 0.10 matches `experiments/004_old_dogs/results.json` (0.102).

## How this compares with published work

`papers/literature-check.md` records that after a treatment stops, the benefit can vanish (Mair et al. 2003, flies), part of it can last (Hahn et al. 2019, mice), or risk can rebound above the untreated level (McCracken et al. 2020, flies: 3.7 times higher hazard, gone after 8 days). Those are animal experiments with thousands of animals and none is in dogs. None was opened in this run. Our rebound sizes are not taken from them; a fly's 8 days has no agreed equivalent in a dog.

The rule from human medicine that after-treatment effects should be estimated from the after-treatment period (Nayak et al. 2018, as recorded in the literature check) matches what we find: the test over the whole trial blends the two periods and can cancel or reverse.

The survivor-mix warning (Bartlett et al. 2020, as recorded in the literature check) is correct in principle and we reproduce its direction: the ratio is above 1 when the treatment truly does nothing after stopping. In this design it is about 1.01 to 1.02, too small to see. That is a statement about this design, not about the general warning.

The literature check did not find, in the searches run, published error rates for a rebound or a partly lasting benefit in an aging-drug withdrawal trial, so there is no published figure to compare these rates with.

## Limits

1. **Every effect size is an assumption**, as are the six-month length of the short rebound and the hidden-differences spread of 0.30. The spread of 0.10 is a tentative fit from experiment 004.
2. **Best case throughout.** The simulated dogs follow the doubling law exactly, nobody drops out, every owner gives the treatment, and all 580 dogs are enrolled. The trial had enrolled far fewer at the last published count (210 in September 2025). Real detection rates would be lower.
3. **The two-period test reports one number for the whole time after dosing.** It cannot tell a short large rebound from a lasting small one, or a benefit that fades slowly from one that half lasts.
4. **ASSUMPTION: the any-difference test** is a two-sided log-rank test at the 5 percent level. The design paper is recorded as naming a Cox model.
5. **Hidden differences are modelled one way only**: a fixed personal factor with a gamma shape, the same in both groups, and a treatment that helps every dog equally. A treatment that helps frail dogs more would change the survivor mix more.
6. **1,000 trials per row.** Differences of 1 to 2 percentage points between rows should not be read.

## Inputs

- Rapamycin trial: level 0.0066, slope 0.29 per year, entry ages 7 to 8 (half), 8 to 9 (quarter), 9 to 10 (quarter), 580 dogs, 12 months of dosing in 36: Coleman et al. 2025, GeroScience 47:2851, as recorded in `interventions/triad_published_inputs.py` and `papers/literature-check.md`. Not re-opened.
- Hidden-differences spread 0.10: `experiments/004_old_dogs/results.json`. Spread 0.30: ASSUMPTION.
- UK level and slope for the older-group side calculation: `dog_lifespan/results.json`.
- Two-period test: `experiments/006_stop_and_watch/two_period.py`. Slope test: `experiments/001_trial_fingerprint/fingerprint.py`. Log-rank test and entry-age mix: `interventions/triad_published_inputs.py`.

## Files

`rebound.py` (about 1.5 minutes), `results.json`, `PREDICTION.md`.
