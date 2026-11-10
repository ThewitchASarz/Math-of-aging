# The research loop

This project advances one question at a time. Every pass through the loop follows the same eight steps, and every pass is logged in `NOTEBOOK.md`, including the ones that fail.

## The eight steps

1. **Pick a question** from the top of `QUESTIONS.md`. Take the highest-ranked question that can be answered with the data and tools available now.
2. **State a hypothesis that could be wrong**, and write down what result would refute it, before running anything.
3. **Do the work in code.** Every number must come from a script in this repository. Fixed random seeds. Results go to a `results.json`, figures to `figures/`, and the series behind each figure to `data/`.
4. **Check the sources.** Every outside number is traced to a page that was actually opened, with its URL. Prefer the journal, agency or sponsor over news coverage. Anything that cannot be traced goes in `KNOWN_GAPS.md`.
5. **Get an independent check.** Someone who did not do the work (a second agent, or a person) re-derives the key numbers with their own code, tries to break the conclusion, and reports what they found.
6. **Record the verdict** in `NOTEBOOK.md`: confirmed, refuted, or inconclusive. Record what the checker changed. Negative results are kept.
7. **Update what depends on it**: the papers, `KNOWN_GAPS.md`, the tests, and the question queue (close the question, add the new ones it raised).
8. **Commit and push**, with tests passing.

## Rules

- No claim without code. No number from memory.
- Models are not measurements. A result about a model says what would follow if the model were true. Say which it is.
- Judgement is labelled as judgement.
- Never describe aging as solved, cured or reversed. No treatment has been shown to slow the rate of aging in humans or dogs.
- Describe evidence only. No medical or veterinary advice.
- A question the loop cannot answer with computation is marked "needs experiment" and stays in the queue, with the experiment that would settle it.

## What the loop can and cannot do

It can sharpen predictions, build and test analysis methods, check consistency between data sets, and find errors. It cannot produce the measurements that only a laboratory or a trial can make. Those are listed in `QUESTIONS.md` under "needs experiment", so that anyone able to run them knows exactly what to measure.
