# Experiment 005: prediction

Written before any code for this experiment was run.

## Question

When the worn-out-cell model is solved directly, with production of worn-out cells rising continuously through life, in place of the frozen-production shortcut used in experiment 002, do that experiment's conclusions about dogs still hold?

## Why it matters

Experiment 002 fitted the worn-out-cell model (Karin and others, 2019) to two dog life tables. It did not solve the model itself. It used a shortcut: the death rate at each age was worked out as if production of worn-out cells were frozen at that age's value. The independent check of experiment 002 then found the shortcut is off by up to 20 percent for dogs. That is as large as the differences between the fits. So three things are in doubt: whether speed alone still fails, how well the two-setting fits really do, and the fitted sizes themselves (for example "dogs clear at about half the human rate").

This experiment removes the shortcut and fits again. It also tests a guess about why the two-setting fits missed by 17 to 28 percent: perhaps the shortcut, not the model, was to blame.

A published paper (Raz and others, 2026, *Nature Aging*) has fitted the same model to four UK dog breeds. This project could not open that paper in this run, so no number from it is used or compared here.

## Method, fixed in advance

- **Model.** The same as experiment 002. The load of worn-out cells X changes as: production (eta times age) minus removal (beta times X divided by (kappa plus X)), plus random swings of size eps. A dog dies when X first reaches the lethal load Xc. X starts at zero at birth and cannot go below zero.
- **Direct solution.** The share of dogs still alive at each age is found by solving the model's own equation for the spread of X across dogs (the Fokker-Planck equation) on a grid, with production rising continuously. No random numbers are needed. The death rate in each year of age is minus the logarithm of the share of dogs alive at the start of the year that survive it.
- **Guard on the solver.** Before fitting, the direct solution is compared with a brute-force simulation of 40,000 dogs (the method of `experiments/002_dog_cell_model/sde_check.py`) at experiment 002's fit-B settings for the US table. The solver is accepted only if the two agree within simulation noise.
- **Human reference settings** (published, Karin and others 2019, as recorded in experiment 002): eta 0.00135, beta 0.15, kappa 0.5, eps 0.142, Xc 17.
- **Data.** The UK dog life table, ages 2 to 17 (Teng and others 2022), and the US table, ages 2 to 16 (Montoya and others 2023), as transcribed in `dog_lifespan/data/`. The same as experiment 002.
- **Fits**, by least squares on the logarithm of the yearly death rate, the same four as experiment 002:
  - A: production speed only.
  - B: production speed and clearance capacity.
  - C: production speed and random swings.
  - D: production speed and lethal load.
- Half-saturation (kappa) is held at the human value, as before.
- "Error" below means the root-mean-square error of the logarithm of the yearly death rate. An error of 0.2 is a typical miss of about 20 percent.

Experiment 002's values, for comparison: fit A 1.81 (UK) and 1.91 (US); fit B 0.28 and 0.20; fit C 0.23 and 0.17; fit D 0.20 and 0.16. In fit B, clearance was 0.51 to 0.53 of the human value and production 3.4 to 3.9 times.

## Predictions

1. **P1. Speed alone still fails.** With the direct solution, fit A has an error above 1.0 on both tables.
2. **P2. The shortcut was not the cause of the misfit.** For fits B, C and D on both tables (six fits), every error with the direct solution is within 0.08 of experiment 002's value, and every one stays above 0.12.
3. **P3. The fitted sizes move a little, and "dogs clear less" survives.** In fit B the clearance capacity is between 0.40 and 0.65 of the human value on both tables, and the production speed is between 2.8 and 4.6 times the human value.
4. **P4. Death records still cannot choose.** On each table the errors of fits B, C and D with the direct solution all lie within 0.10 of each other, and the order of the three from best to worst is the same as in experiment 002 (D best, then C, then B) on both tables.

## What would change our mind

If P1 fails, dogs could be "fast humans" in this model after all, and Theorem 10's conclusion would rest on the shortcut. If P2 fails downward (errors well below 0.12), the model fits dogs better than experiment 002 said, and the misfit was our approximation. If P3 fails, the sizes quoted from experiment 002 should be replaced. If P4's ordering changes, that supports what experiment 002's check said: the ranking was never supported.

These are computational experiments. They show what one mathematical model of one mechanism implies when fitted to published death records. They say nothing about what is inside real dogs; only tissue measurements can do that.
