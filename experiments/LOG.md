# Lab log

One entry per experiment, newest first. Each entry says what was predicted, what happened, and what changed as a result.

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

1. **003: a fading effect.** The rapamycin trial treats for one year. Simulate a slope effect that stops when treatment stops, and measure how often the fingerprint test is fooled.
2. **004: does the doubling law hold for old dogs?** Test for levelling-off in the dog tables above age 14, since experiment 001 depends on it.
3. **005: fit the model itself, not the shortcut.** Refit experiment 002 using simulated death rates, to see whether any ranking of the alternatives survives.
