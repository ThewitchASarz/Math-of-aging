# Experiment 005: result

**Solving the worn-out-cell model directly, without the shortcut, leaves experiment 002's conclusions standing.** Speed alone still cannot turn a human into a dog in this model, and fails by even more than before. The two-setting fits improve a little but still miss by 15 to 22 percent. The fitted sizes of the two-setting fits move by less than a tenth. Death records still cannot say which setting differs. All four predictions were confirmed.

This is a statement about one mathematical model of one mechanism, fitted to published death records. It says nothing about what is inside real dogs. Only tissue measurements can do that.

## What was done

Experiment 002 worked out the death rate at each age as if the production of worn-out cells were frozen at that age's value. Its independent check found that shortcut is off by up to 20 percent for dogs. Here the model's own equation is solved on a grid, with production rising continuously, and the same four fits are run again on the same two dog life tables.

"Error" means the root-mean-square error of the logarithm of the yearly death rate. An error of 0.2 is a typical miss of about 20 percent.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | Speed-only fit still fails: error above 1.0 on both tables | 2.17 (UK), 2.16 (US) | Confirmed |
| P2 | Two-setting fits stay within 0.08 of experiment 002's errors and above 0.12 | Changes of -0.066 to -0.009; lowest error 0.146 | Confirmed |
| P3 | In the clearance fit, clearance is 0.40 to 0.65 of human and production 2.8 to 4.6 times | Clearance 0.49 (UK), 0.52 (US); production 3.60, 3.23 | Confirmed |
| P4 | The three two-setting fits lie within 0.10 of each other; order stays D, C, B | Spread 0.034 (UK), 0.021 (US); order D, C, B on both | Confirmed, with a caveat from the check (below) |

## Main table

Ratios are the dog setting divided by the published human setting. The earlier, shortcut-based values of experiment 002 are in brackets.

| Fit | What differs from humans | Error, UK table | Error, US table | Dog settings |
| --- | --- | --- | --- | --- |
| A | Production speed only | 2.17 (1.81) | 2.16 (1.91) | production ×9.1 to ×9.9 |
| B | Production and clearance capacity | 0.22 (0.28) | 0.17 (0.20) | production ×3.2 to ×3.6, clearance ×0.49 to ×0.52 |
| C | Production and random swings | 0.20 (0.23) | 0.15 (0.17) | production ×5.4 to ×6.1, swings ×1.8 to ×1.9 |
| D | Production and lethal load | 0.18 (0.20) | 0.15 (0.16) | production ×4.7 to ×5.3, lethal load ×0.59 to ×0.61 |

## What this means

1. **Theorem 10's conclusion does not rest on the shortcut.** With the model solved directly, the speed-only fit is worse than before, not better. An error of 2.17 is a typical miss by a factor of almost nine.
2. **The shortcut was not the reason the two-setting fits miss.** The errors fall by 0.01 to 0.07, but every fit still misses by 15 to 22 percent over the full range of ages.
3. **The sizes quoted from experiment 002 hold.** In this model a dog needs faster production plus one of: about half the human clearance capacity, a lethal load about 40 percent lower, or random swings about 1.8 times larger. No setting in the two-setting fits moved by more than 8 percent.
4. **Death records still cannot choose among the three.** They differ in error by 0.02 to 0.03, and all three miss the data in the same way (see the check).
5. **The shortcut itself is poor for dogs and good for humans.** At the fitted dog settings the shortcut and the direct solution differ by up to 15 percent (US, fit B) and up to 32 percent (UK, fit B). At the published human settings they differ by at most 3.5 percent between ages 40 and 85.

## Independent check

A separate agent wrote its own solver and its own fitting code before reading ours. It reproduced all eight errors within 0.002 and all fitted settings within 1 percent, and its own brute-force simulation of 200,000 dogs agreed with its solver (overall ratio 1.002). It found no error in the method. It scored all four predictions as met.

What the check changed:

- **A small slip in our human comparison was fixed.** We had used a two-day time step there, which made each simulated year 364 days long. It now uses one-day steps. The largest gap moved from 3.2 to 3.5 percent.
- **The misses have a shape, and it is the same in every fit.** The model is too low at age 2 (by about 40 to 50 percent), too high at ages 5 to 8 (by about 15 to 28 percent) and too low again at ages 13 to 16. This is a mismatch of curve shape, not noise. The likely cause is deaths in young dogs that have nothing to do with aging, which the model does not include. Experiment 002's follow-up fit with a background death rate pointed the same way.
- **Most of the error comes from the youngest ages.** The check refitted on ages 5 to 15 only, and we then repeated it (an extra, added afterwards, not a test of the predictions). The errors fall to 0.04 to 0.06 (UK) and 0.03 to 0.04 (US). The fitted settings stay close: clearance 0.55 to 0.56 of human, production 3.7 to 4.1 times.
- **P4's ordering is not robust.** The check found the order flips on the US table when ages 2 to 5 are left out (the clearance fit becomes best on ages 6 to 15). So "D best" must not be read as "the lethal-load explanation is the better one". This agrees with what experiment 002's check said: no ranking is supported.
- **The brute-force simulation in experiment 002 read low.** It looked for death only at the end of each half-day step, so it missed dogs whose load touched the lethal level and came back within a step. Against the direct solution it reads about 4 to 6 percent low overall, and about 11 percent low at age 2 by the check's count. Our guard simulation here catches those crossings and agrees with the solver within noise. Experiment 002's statement that the shortcut is "off by up to 20 percent" remains about right; the detail of where it is off shifts.
- **Our guard was a weak test on its own** (40,000 dogs; noise of about 10 percent at age 2). The check's 200,000-dog simulation is the stronger evidence that the solver is right.

## Published work on the same question

`papers/literature-check.md` records that Raz and others (2026, *Nature Aging*) fitted this model to four UK dog breeds and concluded that the production rate is the best single predictor of lifespan differences between species. That paper also warns that different combinations of settings fit almost equally well. Our result agrees with the warning: three different second settings fit about equally. We could not open the paper or its supplement in this run, so no number here is compared with theirs. That comparison stays in the queue as experiment 022.

## Limits

- **One input was not re-checked against its source.** The human settings are those recorded in experiment 002 from Karin and others (2019). Neither we nor the check could open that paper in this run.
- **The life tables were not re-opened either.** The check confirmed that both files name their sources and are internally consistent.
- **Half-saturation is held at the human value** (an assumption carried over from experiment 002). "Dogs clear at about half the human rate" is true only inside fit B, with every other setting held at its human value.
- **The model has no deaths unrelated to aging.** That is the probable reason for the miss at young ages.
- **The UK table is built from deaths in one period, and its youngest rates are rounded to three decimals** (about 3 percent error at ages 2 to 4, by the check's count).
- The two tables give somewhat different settings, so there is no single "dog value".

## Files

`direct_fit.py` (about 3 minutes on 2 cores), `results.json`, `PREDICTION.md`.
