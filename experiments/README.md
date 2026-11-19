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
| 005 | Do experiment 002's conclusions hold when the worn-out-cell model is solved directly, without the shortcut? | Done. Yes: speed alone fails by more (error 2.17); two-setting fits miss by 15 to 22 percent; clearance about half the human value; no ranking supported. The shortcut is off by up to 32 percent at dog settings, 3.5 at human ones. All 4 predictions confirmed. |
| 006 | Can a two-period analysis read a trial in which treatment stops? | Done. Yes: a fading benefit is reported as lasting in 4 to 6% of trials; a lasting one is detected. Blind to modest slowing. All 4 predictions confirmed. |
| 007 | Are breed and size differences large enough to break the slope test? | Done. Size mix: no, but adjust for it (false alarms 9% unadjusted in short trials). Breed differences: still open. 3 confirmed, 1 refuted. |
| 008 | Are differences in risk between breeds large enough to break the slope test? | Done. No: false alarms stay near 5% even at the widest published spread, with or without a baseline per breed. Widely spread hidden differences still break it. 3 confirmed, 1 refuted. |
| 009 | When a benefit fades slowly or only in part after dosing stops, what does the two-period test report, and does a longer or larger trial read it better? | Done. "Lasting" means "still there on average while watched": a slow fade is reported as lasting in 65% of 580-dog, 3-year trials, a permanent half-lasting benefit in 38% (about 92% with 10 years). 5,000 dogs call a quick fade lasting in about 86%. All sizes are assumptions. 3 confirmed, 1 half and half. |
| 010 | How large and how long must a dog trial be to recognise slowed aging? | Done. About 1,300 dogs treated and followed for 10 years for a 20% slowing; no 4-year trial is enough. 2 confirmed, 2 refuted; first headline corrected by the check. |
| 011 | How many dogs does the 10-year trial need with dropout, bending risk and a mix of sizes? | Done. About 1,300 to 1,600 dogs (80% with 5% dropout a year at 1,300; 73% with 10%); 2,000 gives a margin. False alarms stay at 5 to 6%. Dropout rates are assumptions. 2 confirmed, 1 refuted, 1 half and half. |
| 012 | Does an extra setting for hidden differences between dogs repair the slope test, and what does it cost? | Done. Yes: false alarms for a 50% risk cut fall from 19% to 5.5% at an assumed spread of 0.3, with no loss of power found (88.5% against 88.0%). Tested only for a 10-year trial. 2 confirmed, 1 refuted, 1 half and half. |
| 013 | If breeds differ in how fast risk rises with age, does the slope test still work, and what repairs it? | Done. Yes, nearly: even with the whole breed difference put on the rate, a 50% risk cut is misread in about 7% of 1,300-dog, 10-year trials. A level per breed with one shared rate is enough (about 5 to 6%, power 85 to 89%). A rate for every breed is not needed and gives about 7% on its own account. The share on the rate is an assumption. 3 confirmed, 1 refuted. |
| 014 | When dogs leave a trial for reasons linked to health or treatment, do the analyses report effects that are not there? | Done. The slope test holds (5 to 7% false alarms) unless the share of missed deaths drifts over time in one group (11 to 29%). Sick dogs withdrawn unequally (30 against 15%) create a false risk cut in 69% of trials. Every dropout pattern is an assumption. All 4 predictions confirmed. |
| 015 | What happens to the slope test when the treatment takes a year or two to reach full strength, or some owners do not give it? | Done. A delayed true slowing loses nothing (78 to 81%). A 50% risk cut building up over 2 years is misread as slowed aging in 23% of trials; leaving out the first 2 years repairs it (5 to 7%) but power falls from about 80 to about 52%. 20% of dogs never treated: 67 to 69%; 10% of owners stopping each year: 46%. Every delay and share is an assumption. 3 confirmed, 1 half and half. |
| 016 | With their own published designs, what could the two running dog trials detect? | Done. Rapamycin trial as designed: 30% of placebo dogs die; a year of 20% slowed aging shows in 6% of trials; a fading 50% risk cut in 16%. STAY-like: 83% die; a 20% risk cut is detected in about 90%, a 20% slowing recognised in about 30%. 2 confirmed, 2 half and half. |
| 017 | Does the hidden-differences slope test still work in a 4-year trial, with dropout and a mix of sizes, and when the treatment helps frail dogs more? | Done. Yes for the first two (a 50% risk cut is flagged in about 4% of trials, no power lost). No for the third: frail dogs helped more is flagged as faster aging in 30% of 10-year trials (78% with the plain test); sturdy dogs helped more is flagged as slowed aging in 88% (58% plain), with no dog aging more slowly. All sizes are assumptions. 2 confirmed, 2 half and half. |
| 019 | What does a one-year-dosing trial report if risk rebounds after stopping, or only part of the benefit lasts? | Done. Risk halved then 1.5 times untreated: the main test fires in 28% of 580-dog trials and always points to harm; the two-period test sees it in 76%, a short rebound in 30%. A half-lasting benefit is recognised as lasting in 35%. Survivor mix does not matter here. All sizes are assumptions. All 4 predictions confirmed. |
| 032 | How good must a frailty score measured at entry be to stop a benefit that differs between dogs from being read as slowed aging? | Done. Very good: with sturdy dogs helped more, the plain test with score terms still reports slowed aging in 46% of 10-year trials at score quality 0.5, 32% at 0.7, about 8 to 12% at 0.9 (56% with no score). The score works better as a warning: the difference in benefit is noticed in 88% of trials at quality 0.7. With a very good score the adjusted test leans towards faster aging (10% flagged with no effect). All sizes are assumptions. 2 confirmed, 2 half and half. |
| 034 | If a slope test adjusted for a frailty score gives entry age its own setting, does it stop leaning towards "faster aging", and does it then protect against a benefit that differs between dogs? | Done. The lean goes (about 5% flagged with no effect, against about 11%). The protection is weaker than experiment 032 reported: with sturdy dogs helped more, slowed aging is reported in about 40% of trials at score quality 0.7, 21% at 0.9 and 9% at 1.0; five score groups give about 5% at 1.0 and no gain below. The setting costs 5 to 15 points of power. All sizes are assumptions. All 4 predictions confirmed. |

Experiments 018, 020 to 031, 033, 035 and 036 are still in the queue. Experiments 022 and 023 wait on sources that could not be opened. The queue of next experiments is at the end of `LOG.md`.
