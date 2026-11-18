# Lab log

One entry per experiment, newest first. Each entry says what was predicted, what happened, and what changed as a result.

## 2026-11-18: experiment 009, a benefit that fades slowly or in part

- **Question.** In a trial that doses for one year and then watches, what does the two-period test report when the benefit fades slowly or only part of the way, and does a longer or a larger trial read it better?
- **Predicted.** A benefit that half lasts is recognised as lasting in 25 to 45 percent of 580-dog, 3-year trials, over 85 percent with 10 years and over 95 with 5,000 dogs. A lasting tenth is missed everywhere (under 15, and under 50 percent). A slow fade (halving every 2 years) is reported as lasting in 45 to 70 percent of 3-year trials with a ratio of 0.60 to 0.72, and still in 35 to 70 percent of 10-year trials (ratio 0.78 to 0.92). A quick fade (halving every 6 months) is reported as lasting in under 20 percent of 580-dog trials but 50 to 85 percent of 5,000-dog trials.
- **Happened.** Half lasts: 38, 93 and 99 percent (confirmed). A lasting tenth: 9, 25 and 46 percent (confirmed). Slow fade: 65 percent with ratio 0.65; 59 percent with ratio 0.83 over 10 years (confirmed). Quick fade: 18.5 percent with 580 dogs (confirmed) and 87 percent with 5,000 (on the line, just above the predicted range; half and half). Not predicted: 10 years of watching lets a quick fade wash out (7 percent), and in the published design a quick fade and a small permanent benefit look the same (ratios 0.83 and 0.85).
- **Check.** Reproduced independently with separate code, about 40 numbers, 1,500 to 4,000 trials per case; no errors; one decision in about 3,600 shared trials differed. The check put the 5,000-dog quick-fade figure at about 86 percent (on the line), showed one of our figures (24.8 percent) sits about 2 points above its own (chance, traced), and showed the after-dosing ratio is close to an average of the true effect weighted by when dogs die. Its own quick look at splitting the after-dosing time in two: a quick fade is flagged in 66 percent of 5,000-dog trials and 14 percent of 580-dog trials. Neither we nor the check could open the trial's design paper in this run.
- **Changed.** "Lasting" from the two-period test means "still there on average while we watched", not "permanent". The ratio it gives depends on how long the trial watches. More years help (a real partial benefit is recognised, a quick fade washes out); more dogs alone make the average precise without telling fading from lasting. A trial should say this in advance and plan a look at drift over time. All effect sizes and fading speeds are assumptions.

## 2026-11-17: experiment 005, the cell model solved directly

- **Question.** When the worn-out-cell model is solved directly, with production rising continuously, in place of the shortcut of experiment 002, do that experiment's conclusions about dogs still hold?
- **Predicted.** Speed alone still fails (error above 1.0). The two-setting fits stay within 0.08 of their earlier errors and above 0.12. In the clearance fit, clearance is 0.40 to 0.65 of the human value and production 2.8 to 4.6 times. The three two-setting fits stay within 0.10 of each other, in the same order (lethal load best, then swings, then clearance).
- **Happened.** Speed alone: 2.17 and 2.16, worse than the 1.81 and 1.91 of the shortcut. Two-setting fits: 0.15 to 0.22, down by 0.01 to 0.07. Clearance 0.49 to 0.52 of human, production 3.2 to 3.6 times. Spread between the three fits 0.02 to 0.03, same order. All four confirmed. At the fitted dog settings the shortcut is off by up to 15 to 32 percent; at the human settings by at most 3.5 percent.
- **Check.** Reproduced independently with a separately written solver and a 200,000-dog simulation: all eight errors within 0.002, no errors in the method. The check fixed a small slip (a 364-day year in our human comparison), showed the misses have the same shape in every fit and come mostly from ages 2 to 4, and showed the order of the three fits flips on the US table when young ages are left out. It also showed experiment 002's brute-force simulation read about 4 to 6 percent low. Neither we nor the check could open the paper behind the human settings in this run.
- **Changed.** Experiment 002's conclusions no longer depend on the shortcut. The sizes it quoted stand. The shortcut should not be used for dog-like settings. The model's remaining miss points to deaths in young dogs that are not due to aging, which the model leaves out. No ranking of the three explanations is supported; telling them apart still needs measurements in dog tissue.

## 2026-11-17: experiment 014, dropout linked to health or to treatment

- **Question.** In a 1,300-dog, 10-year trial, when dogs leave for reasons linked to their health or to the treatment, do the standard analyses report slowed aging or a cut in risk that is not there?
- **Predicted.** Dropout linked only to treatment, or to health equally in both groups, leaves the slope test at 3 to 8 percent false alarms. Sick dogs withdrawn before death more often in the treated group (30 against 15 percent) make a useless treatment look like a risk cut in over half of trials (risk ratio 0.78 to 0.88) but leave the slope test alone. The slope test is fooled only when the share of missed deaths drifts over time in one group: 6 to 15 percent for a drift to 15 percent, 15 to 40 for a drift to 30. A true 20 percent slowing is still recognised in 65 to 82 percent.
- **Happened.** Slope test 5 to 7 percent in the fair-looking cases. Unequal withdrawal: a false risk cut in 69 percent of trials, middle risk ratio 0.83, slope test 6 percent. Drift: 11 and 29 percent, always pointing to slowed aging. Frail treated dogs leaving: 6 percent. True slowing recognised in 68 to 75 percent (82 with fair dropout). All four confirmed.
- **Check.** Reproduced independently with separate code, 25 numbers, 600 to 3,000 trials per case; no errors. The check showed two of our first figures (8.3 and 8.9 percent false alarms, from 800 trials) were noise: larger runs by both sides give about 6. It also showed the slope test fires in 5 to 7 percent even with fair dropout. No outside source was opened in this run.
- **Changed.** The danger from unfair dropout is mostly to the ordinary "did it cut deaths" result, not to the slope test. A trial should record what happens to every dog that leaves, in both groups, and compare withdrawals between the groups year by year. Every dropout pattern here is an assumption; no published pattern for a dog trial was found.

## 2026-11-17: experiment 019, risk that rebounds after stopping

- **Question.** In a trial that doses for one year and watches for three, what do the standard analyses report when risk rebounds after dosing stops, or when only part of the benefit lasts, and does the mix of surviving dogs bias the after-dosing comparison?
- **Predicted.** With risk halved in the dosing year and then 1.5 times the untreated level, the main test fires in 20 to 45 percent of trials and nearly always points to harm. The two-period test sees that lasting rebound in 60 to 85 percent and a short one (double for six months) in under 35. When half of a halving lasts, the main test fires in 50 to 75 percent and the benefit is recognised as lasting in 20 to 45. Survivor mix does not matter: false alarms under 8 percent and a survivor risk ratio under 1.03.
- **Happened.** Main test 28 percent, every time pointing to harm. Two-period test 76 percent for the lasting rebound, 30 for the short one. Half-lasting benefit: 63 and 35 percent. Survivor mix: 5.8 percent false alarms and a ratio of 1.013 to 1.015. All four confirmed. Not predicted: under a short rebound the main test sits at chance (5.5 percent), because the benefit and the rebound cancel.
- **Check.** Reproduced independently with separate code and 8,000 trials per case, about 30 numbers; no errors. The check showed our survivor ratio and our brute-force test of it measured two different things (same-age dogs against whole groups); both are now computed and reported. Neither we nor the check could open the trial's design paper or the three animal papers in this run.
- **Changed.** A trial that stops dosing should plan the two-period comparison in advance: the test over the whole trial can cancel a benefit against a rebound or read a helpful-while-given treatment as harmful. Survivor mix is cleared as a practical threat for this design. New limits: one after-dosing number cannot tell a short large rebound from a lasting small one, and "no lasting benefit found" is weak evidence at this size. All effect sizes are assumptions; no rebound is known in dogs.

## 2026-11-17: experiment 016, the two real trials with their real designs

- **Question.** When the two dog trials now running are simulated with their own published designs, what could each one detect, and what share of placebo dogs dies?
- **Predicted.** Rapamycin trial as published: a year of 20 percent slowed aging is invisible (under 10 percent); a 50 percent risk cut that fades after the dosing year is detected in 25 to 50 percent and misread by the slope test in over 20. STAY-like trial: 60 to 85 percent of placebo dogs die, a 20 percent risk cut is detected in over 85 percent, a 20 percent slowing recognised in under 35. The rapamycin trial's dogs treated for 10 years: 80 to 95 percent with 1,300 dogs, under 65 with 580.
- **Happened.** Slowed aging in the dosing year: 6 percent (confirmed). Fading 50 percent cut: detected in 16 percent (refuted, too hopeful) and misread in 38 (confirmed). STAY-like: 83 percent of placebo dogs die, the risk cut is detected in 89 to 92 percent, the slowing recognised in 30 to 35 (confirmed, the last on the line). Ten years: 78 percent with 1,300 dogs (refuted, narrowly), 46 with 580 (confirmed). About 30 percent of placebo dogs die in the rapamycin design.
- **Check.** Reproduced independently with separate code, more than 30 numbers; no errors. The check showed the misreading is not an artefact of search limits, and that two trial facts (a second early look after 400 dogs; STAY dosing throughout) are weakly recorded. Neither we nor the check could open the rapamycin trial's design paper in this run.
- **Changed.** The earlier findings hold with the real designs. Two things are new: a benefit that ends with dosing is almost undetectable in the rapamycin design, and the 1,300-dog figure of experiment 010 is the low end, since it falls just short with the rapamycin trial's published inputs. The slope test must not be applied to a trial that stops dosing.

## 2026-11-16: experiment 012, a slope test that allows for hidden differences

- **Question.** If the slope test is given an extra setting for hidden differences between dogs, does it stop misreading a plain cut in risk as a change in the rate of aging, and how much power does it keep?
- **Predicted.** With widely spread hidden differences (0.3) the plain test gives over 15 percent false alarms and the new test 3 to 8; the new test stays at 3 to 8 percent wherever its assumed shape is right; the repair costs at least 5 points of power when it is not needed; with a different shape the plain test is above 15 percent and the new test below 10.
- **Happened.** Plain 19 percent, new 5.5 (confirmed). New test 3.5 to 6.9 percent in all six cases (confirmed). No cost in power: 88.5 against 88.0 percent with no hidden differences, and 81 to 85 percent with them (refuted in its main part, in the useful direction). With a different shape the plain test gives 12.5 percent, not over 15, and the new test 4.5 (half and half).
- **Check.** Reproduced independently with separate code in six cases; no errors; inputs confirmed. The check stopped us calling the higher power in the wide-spread worlds a clean gain, since the plain test is the wrong model there, and showed false alarms can only be stated as "about 5 percent".
- **Changed.** The threat left open by experiments 004 and 008 has a repair for a long trial: add one setting for hidden differences. The plain test's false alarms pointed to faster aging, not slower. Untested: short trials, dropout and sizes with the new test, and treatments that help frail dogs more.

## 2026-11-11: experiment 011, a realistic version of the settling trial

- **Question.** How many dogs does the 10-year trial need once dogs drop out, risk bends at old ages, and the dogs are a mix of sizes?
- **Predicted.** Size mix and bending cost little (80 to 95 percent with 1,300 dogs); 10 percent dropout a year costs more (60 to 78 percent); with 5 percent dropout 1,300 dogs reach 70 to 82 percent and 1,600 are enough; with 10 percent it takes 2,000; false alarms stay at 3 to 9 percent.
- **Happened.** 89 percent with size mix and bending (confirmed). 79 percent with 10 percent dropout (refuted, narrowly: dropout costs less than predicted). With bending and 5 percent dropout, 1,300 dogs reach about 80 percent and 1,600 reach 89 (confirmed). With 10 percent dropout, 1,300 reach about 73 percent and 1,600 sit on the 80 percent line, so 2,000 was too gloomy (refuted, on the edge). False alarms 5 to 6 percent (confirmed).
- **Check.** Reproduced independently with separate code in 14 cases; no errors. The check showed both threshold cells sit within noise of the 80 percent line, so the answer is stated as a range. It could not re-open the paper behind the size table.
- **Changed.** The design of experiment 010 stands but has no margin: about 1,300 to 1,600 dogs, or about 2,000 for a margin against heavy dropout. Dropout rates are assumptions; dropout linked to health is untested and goes to the queue.

## 2026-11-11: experiment 008, differences between breeds

- **Question.** Are differences in risk between breeds large enough to make the slope test misread a risk cut as slowed aging, and does a baseline per breed fix it?
- **Predicted.** Safe at spread 0.15; broken at 0.29 (over 10 percent false alarms) and 0.71 (over 15 percent); a baseline per breed gives 3 to 8 percent and keeps power above 80 percent.
- **Happened.** Safe at every spread: a 50 percent risk cut is misread in about 5 percent of 1,300-dog, 10-year trials even at 0.71 (refuted). A baseline per breed gives 4 to 6 percent and 88 to 91 percent power (both confirmed). A control with widely spread hidden differences of 0.3 still breaks the test (22 percent).
- **Check.** Reproduced independently with separate code; no errors; all 19 rows of the breed table matched to the source paper. The check supplied the reason: five rare, short-lived breeds carry about 82 percent of the spread but are about 3 percent of trial dogs.
- **Changed.** Published breed differences in level of risk are no longer an open threat to the slope test. The "spread" number is a poor yardstick. The open threats are differences inside breeds and breeds that age at different rates; both go to the queue.

## 2026-11-10: experiment 010, the trial that would settle it

- **Question.** How large and how long must a dog trial be to recognise slowed aging?
- **Predicted.** 1,000 dogs for 10 years for a 20 percent slowing; 1,300 dogs for 6 years would reach 55 to 80 percent; 600 dogs for 10 years for a 30 percent slowing; nothing at 4 years.
- **Happened.** 1,000 dogs reach 78 to 79 percent, so the answer is about 1,300 (refuted, narrowly). 1,300 dogs for 6 years reach about half (refuted, on the edge). 600 dogs suffice for 30 percent (confirmed). At 4 years even 5,000 dogs reach only two in three (confirmed).
- **Check.** Ten cells reproduced independently. The check showed our first headline, 1,300 dogs for 8 years, was a lucky draw (77 percent with 1,500 trials) and it was withdrawn.
- **Changed.** The project now has a concrete design to recommend: about 1,300 dogs, treated throughout, followed for 10 years, analysed with a baseline per size class.

## 2026-11-10: experiment 007, differences between breeds and sizes

- **Question.** Are real differences between dogs big enough to break the slope test, and does modelling size fix it?
- **Predicted.** Breed spread above 0.2; ignoring size breaks the test in long trials; a baseline per size class fixes it and keeps power.
- **Happened.** Breed spread 0.71 as published, but 0.29 without one unreliable breed figure (confirmed with a caveat). Ignoring size does not break long trials (refuted) but inflates false alarms to 9 percent in 4-year trials. Size-aware analysis: 4 to 8 percent, 94 percent power (both confirmed).
- **Check.** Reproduced independently. The check traced the French Bulldog problem to the source paper and confirmed the 9 percent with 2,000 trials.
- **Changed.** Recommendation: adjust for size class. Whether breed differences break the test stays open and moves to the queue.

## 2026-11-10: experiment 006, an analysis for stop-and-watch trials

- **Question.** Can a two-period analysis read a trial in which treatment stops?
- **Predicted.** Not fooled by a fading benefit; detects a lasting one; blind to one year of modest slowing.
- **Happened.** Fading benefit reported as lasting in 4 to 6 percent of trials; a persisting 30 percent cut in 65 percent; a year of stopped aging in 56 percent (94 percent with 10 years); 20 percent slowing in 8 percent. All four confirmed.
- **Check.** Reproduced independently; no errors. The check added limits: the test cannot tell slowed aging from a persisting risk cut, a partial fade was not tested, and strong hidden differences between dogs inflate false alarms.
- **Changed.** The misreading found in experiment 003 has a fix. Two limits added to `KNOWN_GAPS.md`.

## 2026-11-10: experiment 004, old dogs and hidden frailty

- **Question.** Does the doubling law hold for old dogs, and does hidden frailty break the fingerprint test?
- **Predicted.** Levelling-off fits better; the plain law is fine to age 16 and fails by age 18; frailty would push false alarms above 10 percent in a long trial.
- **Happened.** Levelling-off fits far better (confirmed). The plain law is 19 percent too high at 16, one year earlier than predicted (refuted), and 80 percent too high at 18 (confirmed). With frailty of the fitted size, false alarms stay near 6 percent (refuted, in the useful direction).
- **Check.** All numbers reproduced independently; no errors. The check corrected our wording (risk bends, it does not level off) and showed that frailty three times larger would break the test.
- **Changed.** "Follow dogs until death" stands, with a condition: a long trial should model differences between dogs. Two items added to `KNOWN_GAPS.md`.

## 2026-11-10: experiment 003, treatment that stops after one year

- **Question.** What can a trial detect when treatment lasts one year, and can a fading benefit fool the fingerprint test?
- **Predicted.** One year of 20 percent slowed aging is undetectable in a 580-dog, 3-year trial; a year of stopped aging is detected about half the time; a fading risk cut is misread as faster aging.
- **Happened.** 8 percent, 50 percent, and misread in 23 percent of 580-dog trials (93 percent of 5,000-dog trials). All four confirmed.
- **Check.** Reproduced independently, including identical decisions on shared simulated trials; no errors.
- **Changed.** Prediction 2 of the theory paper now says the slope comparison is valid only while treatment continues.

## 2026-11-10: experiment 002, the worn-out-cell model fitted to dogs

- **Question.** Which setting of the model must differ between dogs and humans?
- **Predicted.** Production speed alone fails; adding clearance or noise fits within 15 percent; the data cannot choose between them; dogs clear less.
- **Happened.** Production alone fails by a factor of about six (confirmed). Two-setting fits miss by 17 to 28 percent, not under 15 (refuted). Clearance about half the human value in that fit (confirmed). Which setting differs cannot be told from death records (inconclusive as worded).
- **Check.** All fits reproduced independently. The check showed our frozen-production shortcut is off by up to 20 percent for dogs, and that our predicted cell lifetimes were not valid lab predictions. That claim was withdrawn.
- **Changed.** Prediction 1 of the theory paper now states the size of each alternative. `sde_check.py` added. Two items added to `KNOWN_GAPS.md`.

## 2026-11-10: experiment 001, the trial fingerprint

- **Question.** Can the dog trials now running tell slowed aging from delayed death?
- **Predicted.** No: under 30 percent of 580-dog trials and under 50 percent of 1,300-dog trials would recognise a true 20 percent slowing; length matters more than size; false alarms near 5 percent.
- **Happened.** 10 percent and 28 percent; 90 percent with 10 years of follow-up against 33 percent with 5,000 dogs for 3 years; false alarms 3 to 8 percent. All four confirmed.
- **Check.** Reproduced independently with separate code; no errors found; trial sizes and durations confirmed from the sponsors' documents.
- **Changed.** Prediction 2 of the theory paper proposed comparing slopes in the running trials. It now says the comparison needs follow-up until most dogs have died.

## 2026-11-10: experiments 001 and 002 opened

Predictions recorded in `001_trial_fingerprint/PREDICTION.md` and `002_dog_cell_model/PREDICTION.md` and committed before any experiment code was run.

## Next in the queue

1. **022: compare our fitted cell-model settings with Raz et al. 2026.** Raised by the literature check. Open that paper's supplement and compare breed by breed with experiment 002. Waiting on the source: the paper could not be opened on 17 November 2026 (web pages cannot be opened in an unattended run). Take this when the paper can be opened; otherwise go to the next item.
2. **023: between-breed spread from a better table.** Raised by the literature check. Compute the spread of risk between breeds from the 155-breed table of McMillan et al. 2024. Waiting on the source, for the same reason. Take this when the table can be opened; otherwise go to the next item.
3. **013: breeds that age at different rates.** Repeat experiment 008 with breed-level differences in the rate of aging, not only in the level of risk.
4. **015: a delayed or partial effect.** Repeat experiment 011 with an effect that builds up over the first year or two, and with some owners not giving the treatment.
5. **017: the hidden-differences test in harder conditions.** Repeat experiment 012 in a 4-year trial, with dropout and a mix of sizes, and with a treatment that helps frail dogs more than sturdy ones.
6. **018: how the settling trial's size depends on the inputs.** Experiment 016 found 1,300 dogs reach 78 percent with the rapamycin trial's published inputs against 88 with ours. Vary the slope (0.29 to 0.36) and the entry ages separately to find which one matters and how many dogs each case needs.
7. **020: how many dogs a stop-and-watch trial needs.** Experiment 019 found a 580-dog trial sees a 1.25-fold rebound in about a quarter of trials and a half-lasting benefit in about a third. Find the number of dogs and years of watching that bring each to 80 percent, with the rapamycin trial's published inputs.
8. **021: splitting the time after dosing.** Give the two-period test a third period (the first six months after stopping) and measure whether it can tell a short large rebound from a lasting small one, and what that costs in false alarms.
9. **024: how much tracing of withdrawn dogs is enough.** Experiment 014 found that sick dogs withdrawn unequally create a false risk cut, and a drifting share creates false slowed aging. Repeat with the trial finding the date of death for half, or nine in ten, of the withdrawn dogs, and measure how much of each false result remains. Also test whether a year-by-year comparison of withdrawals between the groups would flag the problem.
10. **025: the cell model with deaths unrelated to aging.** Experiment 005 found the directly solved model misses dog death rates mostly at ages 2 to 4, in the same way in every fit. Add a constant background death rate as a fitted setting, solve directly, and measure whether the miss disappears, whether the fitted sizes (clearance about half, lethal load about 40 percent lower, swings about 1.8 times) move, and whether any of the three explanations then stands out.
11. **026: does the after-dosing benefit drift?** Experiment 009 found the two-period test cannot tell a fading benefit from a permanent one, and its check found that separate ratios for the first and second year after dosing flag a quick fade in 66 percent of 5,000-dog trials but 14 percent of 580-dog trials. Test a trend in the after-dosing ratio properly: false alarms, and the dogs and years needed to tell a fade from a permanent benefit of the same average size. Could be run together with 021.
12. **027: partial and slow fades in harder conditions.** Repeat experiment 009 with 5 to 10 percent dropout a year, hidden differences between dogs, risk that bends at old ages (for the 10-year design) and the 210 dogs enrolled at the last published count.
