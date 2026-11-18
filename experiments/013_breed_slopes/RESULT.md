# Experiment 013: result

**Breed differences in the rate of aging, of the size the published life expectancies allow, barely disturb the slope test. Even when the whole difference between breeds is put on the rate (the extreme case), a plain 50 percent cut in risk was misread as a change in the rate of aging in about 7 percent of 1,300-dog, 10-year trials, against the 5 percent the test is designed to have. A level per breed with one shared rate, the simple analysis of experiment 008, brings that to about 5 to 6 percent and has the best power. The heavier repair we planned, a level and a rate for every common breed, is not needed and is slightly worse: it gives about 6.5 to 7 percent false alarms on its own account, even when breeds do not differ at all.** The main prediction (P2) was refuted: we expected 8 to 20 percent. The other three were confirmed, one with a caveat.

This is a computational experiment. It shows what the mathematics and the published breed figures imply under stated assumptions. It does not replace a real trial, and nothing here says any treatment works. "Rate of aging" here means the rate at which death risk rises with age.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | With no true effect, each of the three analyses reports a slope change in 3 to 8% of trials, in all four combinations | 4.0 to 7.1% in all 12 cells (check, 3,000 trials each: 4.9 to 6.6%) | Confirmed |
| P2 | Main case, 50% risk cut: ignoring breed misreads it in 8 to 20% of trials, middle rate multiplier 1.02 to 1.10; a level per breed does not repair it (over 8%); half-and-half case under 10% | 7.2% (check 6.7%), multiplier 1.020 (check 1.025); level per breed 6.1% (check 5.5%); half and half 5.5% (check 5.9%) | **Refuted** in its two main parts: the disturbance is smaller than predicted, and a level per breed is enough. The direction (toward faster aging) and the half-and-half part were right |
| P3 | A level and a rate per breed: 3 to 8% in every case with a 50% risk cut | 6.6 to 7.6% (check 6.4 to 7.2%) | Confirmed by the letter, with a caveat: this analysis runs at about 6.5 to 7% even when nothing is wrong, so it is not a clean repair |
| P4 | A level and a rate per breed recognises a true 20% slowing in 75 to 90% of 10-year trials and 15 to 35% of 4-year trials (main case) | 81.5% and 26.4% (check 82.7% and 27.9%) | Confirmed |

The main case is the breed list without the French Bulldog, with the whole breed difference put on the rate (f = 1), in the 1,300-dog, 10-year design.

## Main table

Share of 1,000 simulated trials in which the slope test reports a change. Each figure near 5 percent is uncertain by about 1.5 points either way, and each figure near 80 percent by about 2.5. Differences of 1 or 2 points between single cells are noise.

| Breed list, share on the rate | Design | Truth | Ignoring breed | Level per breed | Level and rate per breed |
| --- | --- | --- | --- | --- | --- |
| Without the French Bulldog, all on the rate | 10 years | No effect | 4.2% | 4.5% | 6.5% |
| | 10 years | 50% risk cut | 7.2% | 6.1% | 7.6% |
| | 10 years | 20% slowing | 76.5% | 85.3% | 81.5% |
| Without the French Bulldog, half and half | 10 years | No effect | 5.1% | 5.2% | 5.9% |
| | 10 years | 50% risk cut | 5.5% | 4.7% | 6.7% |
| | 10 years | 20% slowing | 81.8% | 89.0% | 85.0% |
| All 19 groups, all on the rate | 10 years | No effect | 4.0% | 4.5% | 5.5% |
| | 10 years | 50% risk cut | 7.1% | 5.7% | 6.6% |
| | 10 years | 20% slowing | 75.3% | 85.9% | 82.2% |
| Without the French Bulldog, all on the rate | 4 years | No effect | 4.6% | 4.2% | 7.1% |
| | 4 years | 50% risk cut | 6.8% | 5.8% | 7.5% |
| | 4 years | 20% slowing | 23.6% | 24.0% | 26.4% |

All designs have 1,300 dogs. About 99 percent of untreated dogs die within the 10-year design and 78 percent within the 4-year one.

The typical (median) estimate of the rate multiplier when the truth is a plain 50 percent risk cut (true value 1.00): 1.02 to 1.05 ignoring breed, 0.99 to 1.00 with a level per breed, 0.99 with a level and a rate per breed.

The check's larger runs (3,000 trials per cell, 24,000 trials pooled over all cells with no true change of slope) give the cleanest summary of false alarms: **about 5.5 percent ignoring breed, about 5.3 percent with a level per breed, about 7.0 percent with a level and a rate per breed.**

## How different the breeds are in this experiment

With the whole difference on the rate, breed rates run from 0.27 per year (Jack Russell Terrier, life expectancy 12.7 years) to 0.56 (English Bulldog, 7.4 years); with the French Bulldog included, to 1.05. That looks like a wide range. But among the dogs that actually enter a trial at ages 7 to 12, the rates are close together: an average of 0.315 with a spread (standard deviation) of 0.041, about 13 percent. The four fastest-aging breeds (Chihuahua, American Bulldog, Pug, English Bulldog) are 5.8 percent of dogs at birth and 1.9 percent at entry. The bulk of the dogs (crossbreds, Labradors, Staffordshire Bull Terriers, terriers, spaniels and collies, about 83 percent at entry) have rates between 0.27 and 0.32.

In the half-and-half case the spread at entry is 0.020.

## Why the main prediction failed

We expected rate differences to act like the widely spread hidden differences that broke the test in experiments 004 and 012 (about 20 percent false alarms). They do not, for the same reason found in experiment 008: the large differences sit in a few rare breeds that mostly do not live to enter the trial.

The check measured the disturbance with noise removed, in one simulated trial of 300,000 dogs with a true 50 percent risk cut. Ignoring breed, the estimated rate multiplier settles at 1.025 where the truth is 1.00. With a level per breed it settles at 0.994. In a 1,300-dog trial the chance scatter of that estimate is about 0.08, three times the bias. A bias of that size should lift false alarms from 5 to about 6 percent, which is what both runs found (7.2 percent in ours, 6.7 in the check's).

## What this means

1. **Breed differences in the rate of aging are not a practical threat to the slope test**, across the range the published breed life expectancies allow, with or without the doubtful French Bulldog figure. The warning at the end of experiment 008, that a level per breed "would not repair" rate differences, was too gloomy for differences of this size.
2. **A level per breed with one shared rate is the analysis to plan.** It kept false alarms near 5 percent in every case and had the highest power for a true slowing (85 to 89 percent in the 10-year design, against 75 to 82 ignoring breed). Ignoring breed costs about 9 points of power when breeds differ in rate.
3. **A level and a rate for every common breed is not worth its cost.** It adds about a dozen settings. The check showed that with 1,300 dogs this alone lifts false alarms to about 6.5 to 7 percent, even in a made-up world where every dog has the same level and rate. It also has 3 to 4 points less power than a level per breed in the 10-year design.
4. **No analysis rescues a 4-year trial.** A true 20 percent slowing is recognised in about a quarter of 4-year trials whichever analysis is used, as in experiments 001 and 008.
5. **The threat that remains is hidden differences inside breeds** (experiment 012), which no table of breeds can measure.

## How this compares with published work

Published studies disagree on whether large dogs die younger because their risk starts higher or because it rises faster (Kraus et al. 2013; Yordy et al. 2020; da Silva & Cross 2023, as recorded in `papers/literature-check.md`). This experiment does not settle that. It shows that for the slope test in a randomised trial the answer matters little: with the whole breed difference on the level (experiment 008), half on each, or all on the rate (here), false alarms stay between about 5 and 7 percent. The literature check did not find published figures for how breed differences affect a slope test in a dog trial. That is a statement about a limited search. No breed-by-breed rates were found in a source we could open, which is why the share on the rate is an assumption. The paper most likely to hold them (Raz et al. 2026) could not be opened in this run.

## Independent check

A separate agent wrote its own simulation and all three analyses from a description of the model before reading our code or results. It ran 3,000 trials for each of the 12 cells with its own seeds, 10,000 trials of the main no-effect case, 6,000 trials of a world with no breed differences, and one trial of 300,000 dogs. It then read our code.

- **It found no errors.** The breed levels and rates matched ours to 9 digits. The groups matched (11 breeds stand alone in the 10-year design and 9 in the 4-year one). Every flagged rate agreed within simulation noise. The largest gap was the 4-year risk-cut cell ignoring breed (6.8 percent in ours, 5.3 in the check's), which is within chance for 1,000 trials.
- **It confirmed the inputs.** The level and rate match `dog_lifespan/results.json`. It reported opening the source page for the breed table (Teng et al. 2022) and matching all 19 rows to the paper's Table 3. It noted that the paper's text gives slightly different counts from its table for some breeds (for example 1,614 against 1,620 Jack Russell Terriers); our file follows the table.
- **It checked the search limits.** No true breed rate lies outside what the analyses can reach, and its own looser limits gave the same results. No trial was skipped.

What the check changed:

- **It showed the heavier analysis is over-eager on its own account.** We had read its 5.5 to 7.6 percent as within the predicted band. The check's larger runs put it at about 6.5 to 7 percent, clearly above 5, and the same in a world with no breed differences. P3 now carries that caveat, and conclusion 3 was added.
- **It supplied the explanation** of why the prediction failed (the 300,000-dog trial and the comparison of bias with chance scatter).
- **It asked that single cells not be over-read.** We now quote the pooled false-alarm rates beside the table.
- It added limits, listed below.

## Limits

- **ASSUMPTION: the share of the breed difference that sits on the rate (f).** It is not measured. The whole-difference case is meant as an extreme. If breeds differed in rate by more than their life expectancies allow under this construction (for example a higher rate offset by a lower level), the disturbance would be larger. That was not simulated.
- **ASSUMPTION: every breed follows the doubling law exactly and matches its published life expectancy at birth**, which is treated as known without error. Life expectancy at birth includes deaths of young dogs, which do not follow the doubling law.
- **ASSUMPTION: the breed mix at birth follows the counts in the published table**, which are the dogs in one UK study, not a census. A trial that enrolled many dogs of short-lived breeds on purpose would see more of the difference.
- **No differences inside a breed**, and "crossbred" (about 32 percent of dogs at entry) is treated as one breed with one rate.
- Entry ages are assumptions. No dropout, treatment effect from the day of enrolment, and the doubling law at all ages, which overstates risk from about age 16 (experiment 004).
- The rule that a breed needs 3 percent of the dogs to get its own group was fixed in advance and not varied. The Shih-tzu (2.9 percent) falls just below it.
- Breed is taken as known. A real trial records breed from owners, with errors.
- The 300,000-dog trial, the world with no breed differences and the pooled false-alarm rates are the check's own figures and are not in `results.json`.
- **Sources not opened in this run:** Raz et al. 2026, McMillan et al. 2024 and the rapamycin trial's design paper (Coleman et al. 2025). The fetch tool's permission request was not answered, as the run was unattended. Kraus 2013, Yordy 2020 and da Silva & Cross 2023 were not re-opened; statements about them rest on `papers/literature-check.md`.

## Files

`breed_slopes.py` (about 3 minutes on 2 cores), `results.json`, `PREDICTION.md`.
