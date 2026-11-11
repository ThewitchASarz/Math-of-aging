# Experiments

This folder is the project's working loop. Each experiment takes one open question, makes a prediction, tests it with code, and records what happened, including when the prediction was wrong.

These are computational experiments. They can show what the mathematics and the published data imply, and what a real trial or lab measurement would need in order to settle a question. They cannot replace that trial or measurement.

## The loop

Every experiment follows the same seven steps, in order.

1. **Question.** One question, taken from the predictions in the theory paper or from `KNOWN_GAPS.md`.
2. **Prediction first.** Before any code runs, write what we expect in `PREDICTION.md`, with numbers, and commit it. The git history shows the prediction came before the result.
3. **Method.** Code in the experiment's folder. Fixed random seeds. Inputs come from this repository's fitted values or from a named source.
4. **Result.** Numbers written to `results.json`, figures to `figures/`.
5. **Independent check.** A separate agent, which did not write the code, re-derives the method, re-runs it, recomputes key numbers its own way, and checks each input against its source.
6. **Verdict.** `RESULT.md` scores each prediction as confirmed, refuted or inconclusive, and lists what the check changed.
7. **Record.** One entry in `LOG.md`, a fast test in `tests/`, and any new open questions added to `KNOWN_GAPS.md`.

## Rules

- A refuted prediction is recorded as plainly as a confirmed one.
- No number enters an experiment without a source or a line of code that produces it.
- Inputs that are assumptions are labelled as assumptions in the code and in the result.
- Nothing here is medical or veterinary advice.

## Index

| No. | Question | Status |
| --- | --- | --- |
| 001 | Can the dog trials now running tell a treatment that slows aging from one that only delays death? | Done. No: 10 to 28% of trials would recognise a true 20% slowing; 90% with 10 years of follow-up. All 4 predictions confirmed. |
| 002 | When the worn-out-cell model is fitted to real dog data, which setting must differ from humans? | Done. Speed alone fails sixfold; which other setting differs cannot be told from death records. 2 confirmed, 1 refuted, 1 inconclusive; 1 claim withdrawn after the check. |
| 003 | What can a trial detect when treatment stops after one year? | Done. One year of 20% slowed aging is detected in 8% of 580-dog, 3-year trials; a fading benefit is misread as faster aging in 23%. All 4 predictions confirmed. |
| 004 | Does the doubling law hold for old dogs, and does hidden frailty break the fingerprint test? | Done. The law overstates risk from age 16 (19%) to 19 (95%); frailty of the fitted size does not break the test. 2 confirmed, 2 refuted. |
| 006 | Can a two-period analysis read a trial in which treatment stops? | Done. Yes: a fading benefit is reported as lasting in 4 to 6% of trials; a lasting one is detected. Blind to modest slowing. All 4 predictions confirmed. |
| 007 | Are breed and size differences large enough to break the slope test? | Done. Size mix: no, but adjust for it (false alarms 9% unadjusted in short trials). Breed differences: still open. 3 confirmed, 1 refuted. |
| 008 | Are differences in risk between breeds large enough to break the slope test? | Done. No: false alarms stay near 5% even at the widest published spread, with or without a baseline per breed. Widely spread hidden differences still break it. 3 confirmed, 1 refuted. |
| 010 | How large and how long must a dog trial be to recognise slowed aging? | Done. About 1,300 dogs treated and followed for 10 years for a 20% slowing; no 4-year trial is enough. 2 confirmed, 2 refuted; first headline corrected by the check. |
| 011 | How many dogs does the 10-year trial need with dropout, bending risk and a mix of sizes? | Done. About 1,300 to 1,600 dogs (80% with 5% dropout a year at 1,300; 73% with 10%); 2,000 gives a margin. False alarms stay at 5 to 6%. Dropout rates are assumptions. 2 confirmed, 1 refuted, 1 half and half. |

Experiments 005 and 009 are still in the queue. The queue of next experiments is at the end of `LOG.md`.
