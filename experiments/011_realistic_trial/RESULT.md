# Experiment 011: result

**The 1,300-dog, 10-year design of experiment 010 holds up better than expected, but it has no margin. With a mix of sizes, risk that bends at old ages and 5 percent of dogs dropping out each year, 1,300 dogs recognise a 20 percent slowing of aging in about 80 percent of simulated trials. With 10 percent dropout a year they reach about 73 percent, and it takes about 1,600 dogs to get back to 80 percent. A trial of 2,000 dogs reaches 88 percent or more in both cases.** The test stays honest throughout: it reports a slope change about 5 to 6 percent of the time when there is none. Two predictions were confirmed, one was refuted narrowly, and one was half confirmed and half refuted.

This is a design calculation under stated assumptions. It says how big a trial would have to be, not that any treatment works. The dropout rates are assumptions, not measurements.

## Predictions scored

| | Prediction | Result | Verdict |
| --- | --- | --- | --- |
| P1 | Size mix and bending, no dropout: 1,300 dogs reach 80 to 95% | 89% (check: 88%) | Confirmed |
| P2 | Size mix and 10% dropout a year: 1,300 dogs reach 60 to 78% | 79% (check: 79%) | **Refuted**, narrowly: dropout costs less than predicted |
| P3a | Realistic scenario: 1,300 dogs reach 70 to 82%, and 1,600 are enough | 79% (check: 82%); 1,600 dogs reach 89% (check: 86%) | Confirmed, on the edge of the band |
| P3b | Harsh scenario: the smallest sufficient trial on the grid has 2,000 dogs | 1,600 dogs reach 80% (check: 81%) | **Refuted**, on the edge: 1,600 sits on the line |
| P4 | In both scenarios, a slope change is reported in 3 to 9% of trials with no effect or a plain 50% risk cut | 5 to 6% in all four cases | Confirmed |

## Main tables

Share of 800 simulated trials in which the slope test recognises a true 20 percent slowing of aging. All trials last 10 years; dogs enter at ages 7 to 12. Each figure is uncertain by about 1 to 1.5 percentage points.

**What each difficulty costs, with 1,300 dogs**

| Scenario | Slowing recognised | Deaths seen | Dogs that drop out |
| --- | --- | --- | --- |
| Size mix only | 94% | 92% | 0% |
| Size mix and bending | 89% | 90% | 0% |
| Size mix and 5% dropout a year | 88% | 73% | 22% |
| Size mix and 10% dropout a year | 79% | 59% | 38% |
| "Realistic": bending and 5% dropout | 79% | 72% | 22% |
| "Harsh": bending and 10% dropout | 75% | 58% | 39% |

**How many dogs are needed**

| Dogs | Realistic | Harsh |
| --- | --- | --- |
| 1,300 | 79% (check: 82%) | 75% (check: 72%) |
| 1,600 | 89% (check: 86%) | 80% (check: 81%) |
| 2,000 | 95% | 88% (check: 88%) |
| 2,600 | 97% | 95% (check: 96%) |

**False alarms, with 1,300 dogs** (share of trials reporting a slope change when there is none)

| Scenario | No effect | Plain 50% risk cut |
| --- | --- | --- |
| Realistic | 6% | 5% |
| Harsh | 6% | 6% |

## What this means

1. **The design of experiment 010 survives, without a safety margin.** Under the realistic scenario 1,300 dogs sit at about 80 percent (79 here, 82 in the check). Under the harsh one they fall to about 73 percent.
2. **Read the answer as a range: about 1,300 to 1,600 dogs.** Around 2,000 dogs gives a margin: 88 percent even with 10 percent dropout a year.
3. **Dropout costs less than the share of dogs lost.** With 10 percent a year, 38 percent of dogs leave before their death is seen, yet power falls only from 94 to 79 percent. A dog that leaves after six years has already contributed six years of information.
4. **Bending and size mix cost a few points each, not more.** A hidden risk factor of the fitted size takes power from 94 to 89 percent.
5. **None of the three makes the test dishonest.** False alarms stay at 5 to 6 percent, and the estimated slowing stays centred on the true 20 percent.

For a real trial this means: plan for dropout in the sample size, and keep following dogs that stop treatment or move away wherever possible, because every death seen counts.

## Independent check

A separate agent wrote its own simulation and its own fitting code from the description, before reading ours. It ran 14 cases with 500 to 2,000 trials each. Every case agreed with ours within simulation noise; the largest gaps were about 3 points in either direction. It confirmed the dropout shares to three digits, verified the formula for the hidden risk factor among surviving dogs by simulating 4 million dogs from age 2, and found no errors in the simulation, the likelihood, the optimiser or the bookkeeping. It confirmed that the script reads its inputs from the repository files named below.

What the check changed:

- **The headline is stated as a range.** Our first reading was "1,600 dogs are enough in both scenarios". The check showed that both threshold cells sit within about one standard error of the 80 percent line (harsh, 1,600 dogs: 80.4 percent here, 80.8 percent in the check with 2,000 trials; realistic, 1,300 dogs: 79 against 82). So we say "about 1,300 to 1,600".
- **False alarms are stated as "5 to 6 percent" and no finer.** The check drew 6 to 8 percent in four cases, then ran a 600-trial calibration that gave 5.3 percent. The differences are noise.
- It added limits 1 and 4 below.

**Not checked:** the checker could not open the paper behind the size table (Montoya et al. 2023), so the five size shares were confirmed against this repository's file, not against the paper again. The file's header records that the paper was opened on 2026-11-10.

## Limits

1. **ASSUMPTION: dropout is unrelated to health and equal in both groups.** Real dropout often is not: owners of sick dogs may leave, or owners may stop a treatment that causes side effects. That could bias the test, not only weaken it. This was not simulated.
2. **ASSUMPTION: dropout rates of 0, 5 and 10 percent a year.** No published rate was used. We tried to open the design paper of the rapamycin dog trial for its planned figure and could not. The literature check in `papers/literature-check.md`, done separately, records other figures from that paper but no dropout rate. The names "realistic" and "harsh" are labels.
3. **ASSUMPTION: the hidden risk factor of experiment 004 applies inside every size class.** Part of the bend in the UK table comes from the size mix itself, so bending is probably counted partly twice. That makes these scenarios slightly too gloomy.
4. **The analysis has the right shape for the simulated dogs** (a doubling law per size group, an effect that is constant from enrolment). A real analysis will fit less well, so these figures are likely on the high side.
5. **ASSUMPTIONS: entry ages 7 to 12; the size mix of trial dogs equals the US population mix; size-class rates of aging are a calibration to life expectancies, not fits to full life tables.**
6. Perfect adherence, and an effect that acts fully from the day of enrolment, are still assumed. A delayed or partial effect was not simulated.
7. The "size mix only" figure (94 percent) is higher than experiment 010's 87 to 90 percent because the two use different mortality inputs (size-class calibration against the all-dog UK fit). Compare rows inside this experiment, not across experiments.

## Inputs

- Size-class level and rates of aging: `dog_lifespan/data/size_class_gompertz_calibration.csv` (this repository's calibration).
- Size mix: `dog_lifespan/data/montoya2023_size_e0.csv` (Montoya et al. 2023, Frontiers in Veterinary Science 10:1082102).
- Spread of the hidden risk factor, 0.102: `experiments/004_old_dogs/results.json`.
- Likelihood and fitting code: `experiments/007_breed_differences/breeds.py`.

## Files

`realistic_trial.py` (about 4 minutes), `results.json`, `PREDICTION.md`.
