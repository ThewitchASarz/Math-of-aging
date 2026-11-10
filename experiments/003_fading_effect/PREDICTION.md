# Experiment 003: prediction

Written before any code for this experiment was run.

## Question

Experiment 001 assumed a treatment that works from enrolment to the end of the trial. The Dog Aging Project's rapamycin trial treats for one year and then watches for two more. What can a trial detect when treatment stops, and can a short treatment fool the fingerprint test?

## Method, fixed in advance

Same simulation and analysis as experiment 001, with three kinds of treatment that last one year:

- **Slowed aging, benefit kept.** The growth of risk is slowed during the treatment year, then returns to normal. The dog stays "younger" by the amount gained. Slowing by 20 percent for one year makes a dog 0.2 years younger for good; stopping aging for the year makes it 1 year younger.
- **Risk cut that fades.** Risk is cut by a fixed fraction during the treatment year and returns fully to the untreated level afterwards.
- **Risk cut that persists.** For comparison: the cut continues after treatment stops.

Designs: 580 dogs followed 3 years (the rapamycin trial's shape), 580 dogs followed 10 years, and 5,000 dogs followed 3 years. Entry ages 7 to 12. At least 400 simulated trials per case.

## Predictions

1. **P1.** 580 dogs, 3 years: aging slowed by 20 percent for one year, benefit kept, is detected by the usual "did it cut deaths" test in fewer than 10 percent of trials.
2. **P2.** Same design: aging fully stopped for the treatment year, benefit kept, is detected in 40 to 75 percent of trials.
3. **P3.** A 30 percent risk cut that fades fools the fingerprint test. It reports a slope change in more than 8 percent of 580-dog, 3-year trials and more than 25 percent of 5,000-dog, 3-year trials, and in most of those the estimated slope is above 1, which would read as faster aging.
4. **P4.** One year of 20 percent slowed aging, benefit kept, is reported as a slope change in fewer than 10 percent of trials in every design.

## What would change our mind

If P1 is refuted, one year of modestly slowed aging is detectable in a trial of this shape, and the rapamycin trial is better placed than experiment 001 suggested.
