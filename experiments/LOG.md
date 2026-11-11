# Lab log

One entry per experiment, newest first. Each entry says what was predicted, what happened, and what changed as a result.

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

1. **012: a slope test that allows for hidden differences.** Add a hidden-difference term to the analysis and test whether it removes the 22 percent false alarms of the control in experiment 008 without losing power.
2. **014: dropout that depends on health or on treatment.** Repeat experiment 011 with sick dogs more likely to leave, and with treated dogs leaving more often, to see whether the slope test reports slowed aging that is not there.
3. **005: fit the cell model itself, not the shortcut.** Refit experiment 002 by solving the model's equation directly.
4. **009: a partial fade.** Repeat experiment 006 with benefits that fade part of the way.
5. **013: breeds that age at different rates.** Repeat experiment 008 with breed-level differences in the rate of aging, not only in the level of risk.
6. **015: a delayed or partial effect.** Repeat experiment 011 with an effect that builds up over the first year or two, and with some owners not giving the treatment.
