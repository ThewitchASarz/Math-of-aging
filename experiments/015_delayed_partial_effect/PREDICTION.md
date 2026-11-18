# Experiment 015: prediction

Written before any code for this experiment was run.

## Question

In the 10-year trial of experiment 011, what happens to the slope test when the treatment takes a year or two to reach full strength, or when some owners do not give it?

## Why it matters

Experiments 010 and 011 found that about 1,300 to 1,600 dogs followed for 10 years could recognise a 20 percent slowing of the rise of death risk. Both assumed the treatment works fully from the first day and that every owner gives it for the whole trial. Neither is likely in a real trial. A trial that plans a slope test needs to know how many more dogs these two problems cost, and whether either one can make the test report slowed aging when the treatment only cuts risk.

"Rate of aging" here means the rate at which death risk rises with age. Nothing here says any treatment works.

## Method, fixed in advance

**The trial.** The "realistic" scenario of experiment 011, unchanged: dogs of five size classes in the US population mix, each class with its own rate of aging (this repository's calibration); a hidden personal risk factor of the size fitted in experiment 004; entry ages 7 to 12; 10 years of follow-up; 5 percent of dogs lost to follow-up each year for reasons unrelated to health; half the dogs assigned to treatment. 1,300 dogs unless stated.

**True effects.** The rise of risk slowed by 20 percent; or risk cut by 50 percent with the rise unchanged; or none.

**Delay (ASSUMPTION).** The effect builds up in a straight line from nothing at enrolment to full strength after D years, with D = 1 or 2. For slowed aging, the slowing itself grows from 0 to 20 percent. For a risk cut, the logarithm of the risk ratio moves in a straight line from 0 to the logarithm of 0.5.

**Owners who do not give the treatment (ASSUMPTIONS).**

- "Never given": 10, 20 or 30 percent of the dogs assigned to treatment never receive it. Their risk is that of untreated dogs.
- "Stopped": every owner in the treated group has a 10 percent chance each year of stopping for good. After stopping, risk rises at the untreated rate again, from where it stood (the gain already made is kept).

No published figure for either was found in a source we could open. All sizes are assumptions.

**Analysis.** Every dog is analysed in the group it was assigned to, whatever its owner did. The slope test of experiment 011 (a level and a rate for each of four size groups; 5 percent level) is used unchanged. It knows nothing about delay or owners.

**A possible repair.** The same slope test applied only to what happens from 2 years after enrolment: dogs that died or were lost before then are left out, and each remaining dog is treated as entering 2 years older.

**Cases.** At least 600 simulated trials per case, fixed seeds. The combined case (1-year build-up and 20 percent never given) is run with 1,300, 1,600, 2,000 and 2,600 dogs.

## Predictions

Experiment 011 found 79 percent for a 20 percent slowing with no delay and every owner giving the treatment.

1. **P1. A delay in a true slowing costs little.** A 20 percent slowing is recognised in 72 to 82 percent of trials with a 1-year build-up and in 66 to 80 percent with a 2-year build-up.
2. **P2. Owners who do not give the treatment cost a lot.** A 20 percent slowing is recognised in 64 to 76 percent of trials when 10 percent of treated dogs never get it, 52 to 66 percent at 20 percent, and 40 to 56 percent at 30 percent. With 10 percent of owners stopping each year it is 45 to 65 percent. In the combined case, 1,300 dogs give 50 to 65 percent, 1,600 dogs give under 78 percent, and the smallest trial on the grid that reaches 80 percent has 2,000 or 2,600 dogs.
3. **P3. A risk cut that builds up slowly is misread as slowed aging.** A 50 percent risk cut that builds up over 2 years is reported as a change in the rate of aging in 20 to 50 percent of trials, with the middle estimate of the rate multiplier between 0.90 and 0.97 (pointing to slower aging). With a 1-year build-up it is 8 to 25 percent. Owners not giving an immediate 50 percent risk cut do less harm: under 12 percent of trials at both 20 and 30 percent never given.
4. **P4. Leaving out the first 2 years repairs the misreading, at a price.** With the analysis started 2 years after enrolment, the risk cut with a 2-year build-up is reported as a slope change in 3 to 8 percent of trials, and no effect in 3 to 8 percent. A true 20 percent slowing that acts from the first day is then recognised in 45 to 68 percent of trials.

## What would change our mind

If P1 is refuted on the low side, delay matters more than thought and trial sizes must grow for it too. If P2 is refuted on the high side, owners matter less than the simple dilution argument says. If P3 is refuted on the low side, a slow start is not a practical threat to the slope test. If P4 fails on false alarms, a different repair is needed; if it fails on the low side for power, the repair is too costly to plan for.

These are computational experiments. They show what the mathematics and this repository's fitted inputs imply under stated assumptions. They do not replace a real trial.
