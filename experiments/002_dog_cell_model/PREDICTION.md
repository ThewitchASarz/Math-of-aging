# Experiment 002: prediction

Written before any code for this experiment was run.

## Question

Theorem 10 of the theory paper says that, in the worn-out-cell model of aging (Karin and others, 2019), speeding up the production of senescent cells cannot by itself turn human mortality into dog mortality. That was argued from two summary numbers. This experiment fits the model directly to the dog life tables.

Which setting of the model must differ between dogs and humans, and by how much?

## Method, fixed in advance

- Model: the published saturating-removal model, with the hazard computed as the inverse of the mean time for the cell load to reach the lethal threshold at each age's production rate.
- Human reference: the published human settings (production slope 0.00135, clearance capacity 0.15, half-saturation 0.5, noise 0.142, threshold 17).
- Data: the UK dog life table, ages 2 to 17 (Teng and others 2022), and the US table, ages 2 to 16 (Montoya and others 2023), as transcribed in `dog_lifespan/data/`.
- Fits, by least squares on the logarithm of the yearly death rate:
  - A: production slope only.
  - B: production slope and clearance capacity.
  - C: production slope and noise.
  - D: production slope and lethal threshold.
- Half-saturation is held at the human value in every fit, because death records cannot pin down more than about two settings.

## Predictions

1. **P1.** Fit A fails: its typical error is above 30 percent (root-mean-square error of the log death rate above 0.3) on both tables.
2. **P2.** Fits B and C both succeed, with typical error below 15 percent.
3. **P3.** The death records cannot choose between B and C: their errors differ by less than 0.03. Telling them apart needs a measurement in dog tissue.
4. **P4.** In fit B, dogs' clearance capacity is lower than the human value.

## What would change our mind

If fit A succeeds, Theorem 10's conclusion does not survive contact with the full data, and dogs could be "fast humans" in this model after all.
