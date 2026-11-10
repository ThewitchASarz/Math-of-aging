# Data sources

Every file under a `data/` folder is one of two kinds:

1. **Transcribed source tables** (17 files): numbers copied from a published source. Each file starts with a comment giving the source URL. They are listed below.
2. **Figure series** (the other 63 files): numbers computed by this project's code, saved so each chart can be rebuilt. These are original to this project and are under CC BY 4.0.

The "reuse terms" column says what is known. **These terms were not individually verified.** Entries marked "to confirm" should be checked before this repository is made public.

| File | Source | Reuse terms |
| --- | --- | --- |
| `human_mortality/data/ssa_2023_period_life_table.csv` | [US Social Security Administration, period life table 2023](https://www.ssa.gov/oact/STATS/table4c6.html) | US federal government work |
| `dog_lifespan/data/ssa2023_us_human_life_table.csv` | Same as above (copy used by the dog comparison) | US federal government work |
| `human_mortality/data/ny2023_cause_shares_by_age.csv` | [New York State Department of Health, Vital Statistics 2023, Table 34a](https://healthweb-back.health.ny.gov/statistics/vital_statistics/2023/table34a.htm) | State government publication; to confirm |
| `dog_lifespan/data/teng2022_uk_dog_life_table.csv`, `teng2022_breed_e0.csv` | [Teng et al. 2022, Scientific Reports](https://pmc.ncbi.nlm.nih.gov/articles/PMC9050668) | Open-access journal (CC BY expected); to confirm |
| `dog_lifespan/data/montoya2023_us_dog_life_table.csv`, `montoya2023_size_e0.csv`, `size_class_gompertz_calibration.csv` | [Montoya et al. 2023, Frontiers in Veterinary Science](https://www.frontiersin.org/journals/veterinary-science/articles/10.3389/fvets.2023.1082102/pdf) | Open-access journal (CC BY expected); to confirm |
| `dog_lifespan/data/mcmillan2024_size_face.csv` | [McMillan et al. 2024, Scientific Reports](https://www.nature.com/articles/s41598-023-50458-w.pdf) | Open-access journal (CC BY expected); to confirm |
| `dog_lifespan/data/causes_by_age.csv`, `oneill2013_causes_by_age.csv`, `oneill2013_bodyweight_coef.csv` | [O'Neill et al. 2013, The Veterinary Journal (repository copy)](https://rvc-repository.worktribe.com/OutputFile/1409027) | To confirm |
| `dog_lifespan/data/scaling_literature.csv` | [de Magalhaes, Costa and Church 2007](https://datawookie.dev/content/blog/2021/2021-04-12-hagr-database-of-animal-ageing-and-longevity/magalhaes-2007.pdf) and other papers named in the file | Individual published values; to confirm |
| `damage_repair/data/saturating_removal_hazard.csv` | Model output using parameters from [Karin et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6889273) | Computed here; parameters are published values |
| `damage_repair/data/clocks_somatic.csv` | Computed from values in [Cagan et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9021023) | Computed here; inputs are published values |
| `damage_repair/data/cancer_dogs.csv` | Computed from the fitted formula in [Nunney 2024](https://escholarship.org/content/qt8w0923ff/qt8w0923ff.pdf) | Computed here; formula is published |
| `interventions/data/itp_mouse_lifespan.csv` | [NIA Interventions Testing Program](https://www.nia.nih.gov/research/dab/interventions-testing-program-itp/supported-interventions) and the papers named in the file | US federal government programme; individual published values |

The evidence tables in `interventions/data/` (`human_evidence.csv`, `dog_evidence.csv`, `other_animal_evidence.csv`) summarise published trial results in this project's own words, with a source URL in each row.

If a source's terms turn out not to allow redistribution, replace the file with a short script that downloads or re-transcribes it, and note that here.
