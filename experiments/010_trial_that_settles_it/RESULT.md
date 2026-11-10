# Experiment 010: result

**To recognise a treatment that slows aging by 20 percent, a dog trial needs about 1,300 dogs treated and followed for 10 years. No trial of 4 years is enough, even with 5,000 dogs.** A 30 percent slowing needs about 600 dogs for 10 years. Two predictions were confirmed and two refuted, and the independent check corrected our first headline.

This is a design calculation under stated assumptions. It says how big a trial would have to be, not that any treatment works.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | With 10 years, the smallest sufficient trial for a 20% slowing has 1,000 dogs | 1,000 dogs reach 78 to 79%, just short; 1,300 reach 87 to 90% | **Refuted**, narrowly |
| P2 | 1,300 dogs for 6 years recognise a 20% slowing in 55 to 80% of trials | 49% (the check measured 55%) | **Refuted**, on the edge |
| P3 | For a 30% slowing, 600 dogs for 10 years are enough | 89 to 92% | Confirmed |
| P4 | With 4 years, nothing on the grid is enough for a 20% slowing | best is 67 to 69%, with 5,000 dogs | Confirmed |

## Main tables

Share of 300 simulated trials in which the slope test recognises the slowing. Treatment is continuous; dogs enter at ages 7 to 12.

**Aging slowed by 20 percent**

| Dogs | 4 years | 6 years | 8 years | 10 years |
| --- | --- | --- | --- | --- |
| 300 | 8% | 16% | 24% | 38% |
| 600 | 10% | 33% | 46% | 58% |
| 1,000 | 16% | 49% | 69% | 79% |
| 1,300 | 22% | 49% | 81%* | 90% |
| 2,000 | 34% | 69% | 92% | 99% |
| 5,000 | 67% | 98% | 100% | 100% |

**Aging slowed by 30 percent**

| Dogs | 4 years | 6 years | 8 years | 10 years |
| --- | --- | --- | --- | --- |
| 300 | 15% | 30% | 48% | 59% |
| 600 | 25% | 51% | 81%* | 92% |
| 1,000 | 34% | 76% | 95% | 98% |
| 1,300 | 44% | 86% | 99% | 100% |
| 2,000 | 58% | 96% | 100% | 100% |
| 5,000 | 94% | 100% | 100% | 100% |

\* These two cells were lucky draws. The independent check, with 1,500 and 400 trials, measured 77 percent and 80 percent. Treat both as just under the line.

## What this means

1. **A trial that could settle whether a drug slows dog aging by 20 percent: about 1,300 dogs, treated throughout and followed for 10 years.** With 8 years it would take roughly 1,500 to 1,600 dogs (estimated from the check; not run). With 6 years, about 5,000.
2. **Length cannot be bought with size at 4 years.** Even 5,000 dogs reach only about two in three.
3. **The last years matter.** By 8 years about 93 percent of dogs have died and by 10 years about 98 percent, yet power still rises from 77 to 87 percent at 1,300 dogs. The oldest deaths carry the most information about the slope.
4. **Against the trials now running:** 1,300 dogs for 4 years, the size of the largest one, gives 22 percent on this grid. The dogs are already enrolled; what is missing is follow-up to death and treatment throughout.

Add the lessons of the earlier experiments: analyse with a baseline per size class (experiment 007), and if treatment stops, use the two-period analysis (experiment 006).

## Independent check

A separate agent wrote its own code before reading ours and reproduced all ten cells it tested. Its values ran slightly lower than ours in most cells. It found no errors in the code, and one in the reading: picking "the smallest size that reaches 80 percent" from 300 trials per cell is fragile, and our first headline of 1,300 dogs for 8 years rested on a lucky draw. With 1,500 trials that cell is 77 percent. The headline above uses its figures.

## Assumptions that matter most

- **Mortality follows the fitted doubling law, and the slowing acts fully from the day of enrolment.** A delayed or partial effect, or mortality that bends at old ages (experiment 004), lowers power.
- **No dogs are lost.** Entry ages 7 to 12, no dropout, no switching between groups and perfect adherence. Real losses mean fewer deaths observed and a larger trial.

## Files

`design_grid.py` (about 75 seconds), `results.json`, `PREDICTION.md`.
