# Experiment 002: result

**Speeding up the production of worn-out cells cannot turn human mortality into dog mortality: that fit misses the dog death rates by a factor of about six.** Letting one more setting differ fixes most of the error. Death records cannot say which setting it is. Two predictions were confirmed, one was refuted and one is inconclusive, and the independent check withdrew one of our claims.

This is a statement about one model of one mechanism. Whether the model describes real dogs is a separate question that only tissue measurements can answer.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | Production-only fit fails (error above 0.3 in log death rate) | 1.81 (UK), 1.91 (US) | Confirmed |
| P2 | Adding clearance or noise brings the error below 0.15 | 0.28 and 0.23 (UK); 0.20 and 0.17 (US) | **Refuted** |
| P3 | Clearance and noise fits differ by less than 0.03 | 0.053 (UK), 0.028 (US) | Inconclusive |
| P4 | In the clearance fit, dogs clear less than humans | 0.51 to 0.53 of the human value | Confirmed |

P2 was too optimistic. Two-setting fits are far better than the production-only fit but still miss by 17 to 28 percent.

## Main table

Error is the root-mean-square error of the logarithm of the yearly death rate; 0.2 means a typical miss of about 20 percent. Ratios are the dog setting divided by the published human setting.

| Fit | What differs from humans | Error, UK table | Error, US table | Dog settings |
| --- | --- | --- | --- | --- |
| A | Production speed only | 1.81 | 1.91 | production ×10 to ×11 |
| B | Production and clearance capacity | 0.28 | 0.20 | production ×3.4 to ×3.9, clearance ×0.51 to ×0.53 |
| C | Production and random swings | 0.23 | 0.17 | production ×5.5 to ×6.3, swings ×1.8 |
| D | Production and lethal load | 0.20 | 0.16 | production ×4.8 to ×5.4, lethal load ×0.60 to ×0.62 |

A follow-up fit, added after seeing these results and so not a test of the predictions, includes a constant background death rate of 0.005 per year (the value fitted to both dog tables). The US errors fall to 0.08, 0.06 and 0.06 for B, C and D.

## What this means

1. **Theorem 10 holds on the full data.** The paper argued from two summary numbers that speed alone cannot produce a dog. Fitting every age of two life tables gives the same answer, and by a wide margin.
2. **The alternatives now have sizes.** In this model a dog needs, besides faster production, one of: about half the human clearance capacity, a lethal load about 40 percent lower, or random swings about 1.8 times larger.
3. **Death records cannot choose among them.** The three fits differ in error by 0.03 to 0.08, and the model's own approximation error is about 0.12 to 0.2 (see below). The ranking of B, C and D is within that error.
4. **A three-setting fit is not identified.** Fit E ran to the edge of its allowed range.

## Independent check

A separate agent wrote its own code for the model and the fitting before reading ours and reproduced every error and ratio in the table. It confirmed the human settings against the source paper. It then found two problems, both now recorded here.

**The shortcut is only approximate for dogs.** We compute the death rate at each age as if production were frozen at that age's value. Simulating the model itself with production rising continuously gives death rates 8 to 18 percent lower than the shortcut at ages 2 to 7 and 10 to 16 percent higher at ages 15 to 16. `sde_check.py` reproduces this (typical error 0.12, largest 0.20). That is as large as the differences between fits B, C and D, so no ranking among them is supported.

**One claim withdrawn.** Our first draft predicted how long worn-out cells last in a young dog under each fit (about 6 to 7 days under B against 3 days under C and D). The formula behind those numbers assumes the cell load sits near its noise-free value. In this model the random swings are much larger than that value, so the numbers are not what a lab would measure. They remain in `results.json` for transparency and should not be used. The qualitative ordering, slower clearance under fit B, is expected to hold but was not checked.

The check also noted that the bootstrap intervals in `results.json` cover sampling error only and are far too narrow to represent the real uncertainty.

## Files

`dog_cell_model.py` (about 70 seconds), `sde_check.py` (about 10 seconds), `results.json`, `sde_check.json`, `PREDICTION.md`.
