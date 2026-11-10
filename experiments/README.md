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

The queue of next experiments is at the end of `LOG.md`.
