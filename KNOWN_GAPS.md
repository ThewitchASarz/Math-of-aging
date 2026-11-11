# Known gaps

Things in this project that are unverified, unresolved or approximate. Read this before relying on any number.

## Not verified against a primary source

- **Dog somatic mutation rate.** The project uses a prediction (200 to 270 substitutions per genome per year) from the inverse-lifespan relationship in Cagan et al. 2022. The literature check (`papers/literature-check.md`) read the measured value in the paper's Supplementary Table 3 as 249.5 per year, inside that range. The project's own code has not been changed to use it, and the table was read by one research agent only.
- **Share of children expected to reach 100 (Olshansky et al. 2024).** Two readings of the paper disagreed on the exact figure (5.1% or 5.3% of girls). The report quotes only the paper's stated ceiling: unlikely to exceed 15% of women and 5% of men this century.
- **TAME metformin trial status.** That it has not enrolled anyone rests on the funder's fundraising page and on news reports from September 2026, not on a trial registry.
- **Loyal LOY-002 approval status.** "No conditional approval as of 10 November 2026" is an absence of evidence in the sources found. The remaining FDA section is described as manufacturing by trade press, not by the company's own page.
- **SPRINT and EMPA-REG OUTCOME hazard ratios** (0.73 and 0.68). Taken from summary pages and not re-checked against the journal papers. SELECT (0.81) and the statin meta-analyses (0.90 in 2010; 0.88 in 2005) were checked.
- **US death rates above age 95.** The Social Security table's rates at these ages are actuarial estimates. Results about the late-life plateau and the oldest attainable age reflect those estimates.
- **Cause-of-death shares by age.** The cause-elimination calculation uses New York State 2023 shares applied to US death rates. A rough check against US national shares (NVSR vol. 74 no. 10) gave the same conclusions, but the calculation was not redone age by age.
- **Mouse single-gene results** (Ames dwarf, Snell dwarf, growth-hormone-receptor knockout, klotho, telomerase). Magnitudes were not verified and are not used in any calculation.
- **Karin et al. 2019 equation form.** The parameter values were confirmed from the paper. The exact form of the equation appears in a figure that could not be read; the form used matches the paper's text and Alon's lecture notes.
- **Frailty-index growth rates** in humans and dogs, and the **historical trend in life expectancy at age 80**: not retrieved.
- **Large-animal safety data for partial reprogramming**: none found.

## Unresolved conflicts

- **Power of the dog rapamycin trial (TRIAD): explained, with one loose end.** The design paper's "greater than 80% power" assumes a risk cut that lasts, with the rate at which risk rises left unchanged, and it is sized for a large effect. Re-running the calculation with the paper's own published inputs (`interventions/triad_published_inputs.py`) gives about 29% power for a 20% risk cut over 3 years, about 80% for a 38% cut, 7% if the cut lasts only for the dosing year, and 9% for a 20% slower rise of risk. Our earlier quoted figure of about 18% was a rough calculation. The code in this repository gives about 23% (`interventions/trial_design.py`: every dog enters at exactly age 7, so 21% of placebo dogs die in 3 years, against 30% with the paper's inputs) and 27 to 43% (experiment 001, depending on entry ages). All of these are far below 80%, and none contradicts the paper, which is sized for a larger effect. Loose end: the paper's two stated effect sizes (6% of mean lifespan; 11% in survival within 3 years) are not explained in the paper, and our reading of them is an interpretation. The trial doses for 12 months of its 36; experiments here that assume treatment throughout do not describe it.
- **Dog-to-human aging speed.** The ratio is 3.3 to 3.7 depending on the fitting method. Statements of "about 3.5" use the middle of that range.

## Calibrated or approximate, not fitted

- **Dog size-class mortality parameters** are calibrated to published life expectancies, one number per size class, because size-specific life tables were not available.
- **The within-dog scaling of lifespan with body mass** (exponent about -0.11) rests on five size-group points and assumed representative masses.
- **"Slowing aging by 10%"** adds 3 to 5 years depending on how slowing is defined.
- **Theorem 11** in the theory paper is an asymptotic approximation (within 2% at ages 40 to 60 for the published human parameters). The other theorems are exact within their stated assumptions.

## Found by the experiments

- **The frozen-production shortcut for the worn-out-cell model is off by up to 20 percent for dogs** (`experiments/002_dog_cell_model/sde_check.py`). Any result that depends on differences smaller than that, including the ranking of which setting differs between dogs and humans, is not supported. The same shortcut is used in `damage_repair/saturating_removal.py`, where its error for humans was measured at about 7 percent in the doubling time.
- **Predicted senescent-cell lifetimes in dogs were withdrawn** as lab predictions; see `experiments/002_dog_cell_model/RESULT.md`.
- **The fingerprint test is valid only while treatment continues** (experiment 003). A benefit that fades after treatment stops is misread as faster aging. No analysis for stop-and-watch trials has been built yet.
- **The doubling law overstates dog death rates from about age 16** (experiment 004): by 19 percent at 16 and by 50 to 95 percent at 17 to 19 when fitted to younger ages. Results here that concern dogs older than 15 lean on it.
- **Hidden differences between dogs.** The fingerprint test holds up under frailty of the size fitted to the UK table, but not under frailty three times larger. Differences between breeds in level of risk were tested in experiment 008 and do not break it.
- **Stop-and-watch trials: use the two-period analysis of experiment 006, and know its limits.** It cannot tell slowed aging from a risk cut that persists; a partial fade was not tested; and "no lasting benefit" is not evidence against modest slowing (power 8 percent for a year of 20 percent slowing in 580 dogs).
- **Ignoring size class inflates the slope test's false-alarm rate to about 9 percent in 4-year trials** (experiment 007). A baseline per size class fixes it.
- **Breed differences in level of risk do not break the slope test** (experiment 008), at any published spread from 0.15 to 0.71. This rests on assumptions: breeds differ only in level, share one rate of aging, and every dog of a breed has the same risk. The French Bulldog figure remains unreliable, but the result does not depend on it.
- **The "spread" number is a poor yardstick for the slope test** (experiment 008). A spread of 0.3 breaks the test when it is spread across all dogs and does nothing when it sits in a few rare, short-lived breeds. Statements in experiments 004 and 007 that compare a spread with "0.3 broke the test" should be read with that in mind.
- **Hidden differences inside breeds are unmeasured.** No source for their size has been found. If they are widely spread with a spread near 0.3, the slope test gives about 22 percent false alarms for a 50 percent risk cut in a 10-year trial. No analysis that allows for them has been built yet.
- **Breeds that age at different rates have not been simulated.** A baseline per breed would not repair that.
- **The trial-size result of experiment 010 was idealised; experiment 011 added dropout, bending risk and a mix of sizes.** The answer becomes a range, about 1,300 to 1,600 dogs, with about 2,000 for a margin. It still assumes perfect adherence and an effect that acts fully from enrolment.
- **Dropout rates in experiment 011 are assumptions (0, 5 and 10 percent a year).** No published dropout rate for a dog trial was found. The experiment could not open the design paper of the rapamycin dog trial (TRIAD, Coleman et al. 2025); `papers/literature-check.md` records other figures from that paper but no dropout rate.
- **Dropout linked to health or to treatment has not been simulated.** Experiment 011 assumes dogs leave for reasons unrelated to health, equally in both groups. If sick dogs or treated dogs leave more often, the slope test could be biased, not only weaker.
- **Experiment 011 applies the hidden risk factor of experiment 004 inside every size class.** Part of the bend in the UK table comes from the size mix, so the bend is probably counted partly twice.
- **The UK dog table is built from deaths in one period**, not by following dogs. Its shape at the oldest ages mixes aging with breed mix and recording errors.

## Judgement, not data

- The ranked estimate of years gained per treatment (`interventions/monte_carlo.py`). Its inputs are the authors' probabilities and ranges, anchored on the evidence tables.
- The bottleneck ranking and scenario outputs in `escape_velocity/bottleneck_voi.py`. These restate the authors' assumptions and should not be read as forecasts.
- Every input marked `ILLUSTRATIVE` in the code.

## Before making this repository public

- **Data reuse terms have not been checked source by source.** `DATA_SOURCES.md` lists every transcribed table and what is known about its terms. Tables from government sources and open-access journals are expected to be reusable with attribution; the others need a check.
- **Novelty: a literature search has been run, an expert review has not.** `papers/literature-check.md` records what was found. Already published: the escape-velocity threshold (de Grey 2004; Missov & Lenart 2011; Goldstein & Wachter 2006), the Gompertz life-expectancy formula and the small gain from removing causes of death (Missov 2013; Keyfitz 1977), frailty bending the mortality curve (Beard 1959; Vaupel, Manton & Stallard 1979), the worn-out-cell model fitted to dog breeds from the same UK data (Raz et al. 2026, Nature Aging), and the difference between lowering the level of risk and slowing its rise, with stop-and-watch designs to separate them (Mair et al. 2003; Simons et al. 2013). Not found in the searches run: how often a dog or human trial would recognise a slower rise of risk, error rates for a fading benefit being misread, and frailty or between-breed variance figures for dogs. "Not found" is a statement about a limited search. No theorem or finding here should be described as original without review by an expert.

## Raised by the literature check

Work suggested by `papers/literature-check.md`. Items marked **experiment** are candidates for the queue in `experiments/LOG.md` and should be taken before older queue items.

- **Experiment: the two real trials with their real designs.** Repeat experiments 001 and 010 with the rapamycin trial's published inputs (level 0.0066, slope 0.29, entry ages 7 to 10, 12 months of dosing in 36 months of follow-up, early looks after 200 and 400 dogs) and with a design like Loyal's STAY study (1,300 dogs aged 10 or older, nearly all sizes, about 4 years, treatment throughout). State the share of placebo dogs that die in each case.
- **Experiment: risk that rebounds after stopping.** Published switch experiments show three outcomes after a treatment stops: the benefit vanishes (Mair et al. 2003), part of it lasts (Hahn et al. 2019), or risk rebounds above the untreated level (McCracken et al. 2020). Experiment 006 tested only the first and a fully lasting benefit. Add the rebound, and the survivor-selection problem described by Bartlett et al. 2020.
- **Experiment: compare our fitted cell-model settings with Raz et al. 2026.** That paper fitted the same model to four UK breeds. Open its supplement and compare breed by breed. Its warning applies here: different settings fit almost equally well.
- **Experiment: between-breed spread from a better table.** Compute the spread of risk between breeds from the 155-breed table of McMillan et al. 2024, which counts dogs still alive, in place of the deaths-only UK table.
- **Wording.** "Rate of aging" in this repository means the Gompertz slope, the rate at which death risk rises with age. Whether that slope measures biological aging is disputed (Zhang & Gems 2026). The papers have not yet been reworded.
- **Credits not yet added inside the two papers.** The theory paper now carries a prior-work note at the top, but individual theorems are not yet cited one by one.
- **The late-life bend and the frailty figure for dogs are tentative.** The UK table counts only deaths, has 78 deaths at ages 19 to 20, and bends at a different age from the US and Japanese tables. The slope itself is better supported: two national tables, the rapamycin trial's own 0.29, and an independent disease analysis (Salt et al. 2023) all fall between about 0.29 and 0.36 per year, and our fit moves by about 0.03 with the choice of ages.
- **The lifespan with risk frozen at its age-30 level** (479 and 1,138 years here) has no published source that we could find; it is this project's arithmetic.
- **Trial facts corrected.** 580 dogs is the rapamycin trial's target; 210 were enrolled at the last published count (September 2025), and the trial is due to end in November 2029. No dog longevity drug is approved as of 10 November 2026 on the company's own pages and the trade press; the FDA's own pages could not be opened.
