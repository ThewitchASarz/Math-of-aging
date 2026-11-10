# Lab log

One entry per experiment, newest first. Each entry says what was predicted, what happened, and what changed as a result.

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

1. **005: fit the model itself, not the shortcut.** Refit experiment 002 by solving the model's equation directly, to see whether any ranking of the alternatives survives.
2. **006: an analysis for stop-and-watch trials.** Build and test a model with separate on-treatment and after-treatment effects, so a fading benefit is not misread as faster aging.
3. **007: differences between breeds.** Measure how much dogs differ by breed and size in the published tables, and test the fingerprint test against differences of that size.
