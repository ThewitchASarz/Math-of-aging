# Question queue

Ranked by how much each answer would matter for extending healthy lifespan, among questions this project can actually work on. The loop takes the top open question that can be answered with data and tools available now. See `PROTOCOL.md` for the steps and `NOTEBOOK.md` for the record.

## Open: can be worked on now

| # | Question | Why it matters | First step |
| --- | --- | --- | --- |
| Q5 | When applied to real, published lifespan data, does the fingerprint test say rapamycin slowed aging or only delayed death in mice? | It is the first use of the method on real data, and the mouse result is the basis for the dog rapamycin trial. | Get the public survival data from the NIA Interventions Testing Program and run `fingerprint/fingerprint.py` on it. |
| Q3 | Why does our dog fit of the senescent-cell model need something beyond faster production, when Raz et al. 2026 (Nature Aging) report that production rate is what best predicts lifespan across species, dogs included? | Our Experiment 2 may be explained, or contradicted, by published work we have not read. | Read the paper's dog parameters and compare them with `dog_model_fit/results.json`. |
| Q4 | Does the fingerprint test stay reliable when the real world is messier: a slowed-aging effect on top of a background hazard, mixed breeds and sizes, dropout, entry spread over years, large dogs only? | A real trial analysis would face all of these. The check so far covered only messy baselines with no slope effect. | Extend the simulator and measure false-alarm and power under each. |
| Q6 | What does removing cancer or heart disease add when the calculation uses US national cause-of-death shares age by age? | Removes one known gap. Expected to change little. | Transcribe NVSR vol. 74 no. 10 and rerun `human_mortality`. |
| Q7 | Do aging clocks pass the marker test on published numbers? A treatment that shifts a clock by d standard deviations should change death risk by the clock's hazard ratio raised to the power d. | It is the cheapest check of whether a clock can stand in for lifespan in trials. | Collect published clock shifts and mortality results from the same trials. |
| Q8 | Can breed-level data turn the dog size-class calibration into a real fit? | Would test directly whether size acts on the rate of aging. | Find a breed-by-breed table of adult weight and age at death. |
| Q9 | Giant dog breeds are reported to lose LINE1 DNA methylation about 35% faster than small breeds (Science, 2026; seen only in news coverage). Is that consistent with their faster rate of aging in our size-class fit? | A measured molecular rate that could be compared with a fitted mortality rate. | Read the paper itself; compare the reported rate ratio with the size-class rates of aging. |

## Open: needs an experiment (cannot be answered here)

| # | Question | The measurement that would settle it |
| --- | --- | --- |
| E1 | In dogs, which differs from humans: clearance of senescent cells, the size of random swings, the lethal load, or the saturation point? | Senescent-cell turnover measured in dog tissue at several ages, as has been done in mice. Experiment 2 gives the model's value for each candidate. |
| E2 | Do the dog drugs now in trials slow aging or delay death? | The survival data from Loyal's STAY study and the Dog Aging Project's TRIAD trial, run through the fingerprint test. |
| E3 | Can any blood or DNA marker stand in for lifespan? | Marker measurements in stored samples from finished trials that counted deaths. |
| E4 | Is any kind of aging damage beyond repair? | Animal studies that remove one kind of damage at a time and measure how much the rise in death risk flattens. |
| E5 | Do effects found in mice carry over to dogs and people at the same hazard ratio? | One drug tested for lifespan in all three species. |

## Closed

| # | Question | Verdict | Notebook entry |
| --- | --- | --- | --- |
| Q1 | Can the dog trials now running tell slowed aging from delayed death? | A trial the size of STAY can, if the effect is large; the rapamycin trial as designed mostly cannot | 2026-11-10, entry 1 |
| Q2 | In the senescent-cell model, is a dog a human with faster damage production? | No, within the model; but the data cannot say what else differs | 2026-11-10, entry 2 |
