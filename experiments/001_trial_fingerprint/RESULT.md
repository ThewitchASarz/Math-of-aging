# Experiment 001: result

**The dog trials now running are too short to tell slowed aging from delayed death.** With 3 to 4 years of follow-up, a treatment that truly slowed aging by 20 percent would be recognised as a slope change in 10 to 28 percent of trials. Following the same dogs for 10 years raises that to 90 percent. All four predictions were confirmed.

These are simulated trials. They say what the planned designs could detect under stated assumptions. They say nothing about whether any drug works.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | 580 dogs, 3 years: a true 20% slowing is recognised in under 30% of trials | 10% | Confirmed |
| P2 | 1,300 dogs, 4 years: recognised in under 50% of trials | 28% | Confirmed |
| P3 | 1,300 dogs for 10 years beats 5,000 dogs for 3 years | 90% against 33% | Confirmed |
| P4 | A pure risk cut is misread as a slope change in 3 to 8% of trials | 3.3 to 7.5% across 20 cases | Confirmed |

## Main table

Share of 400 simulated trials. "Usual test" asks only whether the treatment cut deaths. "Fingerprint" asks whether the slope of risk against age changed.

| Design | Deaths in trial | Usual test detects a 20% slowing | Fingerprint recognises a 20% slowing | Fingerprint recognises a 30% slowing |
| --- | --- | --- | --- | --- |
| 580 dogs, 3 years, entry ages 7 to 12 | about 235 | 12% | 10% | 14% |
| 580 dogs, 3 years, entry ages 7 to 9 | about 155 | 12% | 10% | 9% |
| 1,300 dogs, 4 years, entry ages 10 to 14 | about 1,050 | 50% | 28% | 56% |
| 5,000 dogs, 3 years (hypothetical) | about 2,020 | 62% | 33% | 69% |
| 1,300 dogs, 10 years (hypothetical) | about 1,280 | 99% | 90% | 100% |

For scale, slowing aging by 20 percent from age 7 adds about 0.7 years to a dog's life in this model; a 20 percent risk cut from age 7 adds about 0.5 years.

## What this means

1. **A drug that really slows aging can look like a failure in a short trial.** Slowing the growth of risk changes little at first; the benefit builds up with time. In a 580-dog, 3-year trial a true 20 percent slowing reaches ordinary statistical significance about one time in eight.
2. **Length beats size.** The trial must run until most dogs have died for the slope to be measured. Adding dogs to a short trial helps much less.
3. **How precisely can the slope be measured?** When the true slope multiplier is 0.8, the estimate ranges (5th to 95th percentile) from 0.22 to 1.46 with 580 dogs over 3 years, from 0.59 to 1.04 with 1,300 dogs over 4 years, and from 0.69 to 0.89 with 1,300 dogs over 10 years.
4. **Prediction 2 of the theory paper needs amending.** It proposed comparing slopes in the trials now running. That comparison needs follow-up to death, not the planned 3 to 4 years. The paper has been updated.

## Independent check

A separate agent wrote its own simulation and likelihood code before reading ours and reproduced every number in the table within Monte Carlo noise (for example 0.09, 0.28, 0.35 and 0.90 against our 0.10, 0.28, 0.33 and 0.90). It found no errors in `fingerprint.py`. It confirmed from the sponsors' documents that the rapamycin trial plans 580 dogs aged 7 and over, treated for 1 year and followed for 2 more, and that Loyal's trial has 1,300 dogs aged 10 and over for 4 years.

## Assumptions that matter, and which way they push

- **The effect starts at enrolment and lasts for the whole trial.** The rapamycin trial treats for only one year, so its real power to see a slope change is lower than shown here.
- **Mortality follows the fitted doubling law exactly.** If real dog mortality levels off at the oldest ages, the slope test would be biased. Not simulated.
- **Entry ages are assumed.** Older dogs at entry give more deaths and more power.
- **All-dog UK mortality is used.** The rapamycin trial enrols dogs of 20 kg and over, which die sooner (more power); healthy-volunteer selection works the other way.
- **No dogs drop out.** Optimistic.
- **Baseline mortality is re-estimated in each trial**, as a real analysis would. Treating it as known raises power by 5 to 10 points.

## Files

`fingerprint.py` (about 75 seconds), `results.json`, `PREDICTION.md`.
