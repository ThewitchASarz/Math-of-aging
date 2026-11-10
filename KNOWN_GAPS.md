# Known gaps

Things in this project that are unverified, unresolved or approximate. Read this before relying on any number.

## Not verified against a primary source

- **Dog somatic mutation rate.** The measured value is in the supplementary tables of Cagan et al. 2022, which could not be opened. The project uses a prediction (200 to 270 substitutions per genome per year) from the paper's inverse-lifespan relationship.
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

- **Power of the dog rapamycin trial (TRIAD).** Our calculation gives about 18% power to detect a 20% cut in death risk with 580 dogs over 3 years, if death were the only endpoint. The trial's design paper reports more than 80% power to detect a 6% difference in mean lifespan. The two use different endpoints and assumptions and have not been reconciled.
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
- **Hidden differences between dogs.** The fingerprint test holds up under frailty of the size fitted to the UK table, but not under frailty three times larger. Differences between breeds have not been measured against that threshold.
- **The UK dog table is built from deaths in one period**, not by following dogs. Its shape at the oldest ages mixes aging with breed mix and recording errors.

## Judgement, not data

- The ranked estimate of years gained per treatment (`interventions/monte_carlo.py`). Its inputs are the authors' probabilities and ranges, anchored on the evidence tables.
- The bottleneck ranking and scenario outputs in `escape_velocity/bottleneck_voi.py`. These restate the authors' assumptions and should not be read as forecasts.
- Every input marked `ILLUSTRATIVE` in the code.

## Before making this repository public

- **Data reuse terms have not been checked source by source.** `DATA_SOURCES.md` lists every transcribed table and what is known about its terms. Tables from government sources and open-access journals are expected to be reusable with attribution; the others need a check.
- **Novelty has not been checked.** No theorem should be described as original without a literature review by an expert.
