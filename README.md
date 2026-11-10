# The Math of Aging: Humans and Dogs

A verified mathematical description of how humans and dogs age: mortality models fitted to life tables, models of what goes wrong inside the body, a review of the evidence for treatments, the mathematics of "escape velocity", 19 theorems with proofs, and 5 predictions that an experiment could confirm or kill.

**This is not a cure for aging and does not claim one.** The mathematics says what a treatment for aging would have to do and how to recognise one in trial data. No treatment has yet been shown to slow the rate of aging in humans or dogs.

**Nothing here is medical or veterinary advice.** The treatment sections describe published evidence. They do not recommend what any person or animal should take.

## The key numbers

| | Humans (US, 2023) | Dogs (UK and US) |
| --- | --- | --- |
| Risk of death doubles every | 7 to 9 years | 2.0 to 2.3 years |
| Life expectancy at birth | 75.8 (men), 81.1 (women) | 11.2 (UK), 12.7 (US) |
| Gain from halving every death rate | 7.3 to 8.4 years | about 2 years |
| Pace of progress needed to outrun aging | death rates falling about 8% a year | about 27% a year |

- Removing all cancer deaths adds 2.3 to 2.7 years to US life expectancy; removing all heart-disease deaths adds 2.4 to 3.0.
- If the risk of death stayed at its age-30 level, the average life would run about 480 years (men) and 1,140 years (women).
- US death rates at older ages have fallen 0.5 to 1.0% a year, which is 8 to 17 times slower than the pace needed to outrun aging.
- Everything known today, combined, might add about 5 years (range 2 to 10) to a 50-year-old and about 1 year to a 7-year-old dog. **This estimate is built on judgement calls, not measurements.**

## What is in this repository

| Folder | Contents |
| --- | --- |
| `papers/` | The plain-language report and the theory paper (19 theorems, proofs, predictions) |
| `human_mortality/` | Fits to the US 2023 life table; what moves life expectancy; limits to lifespan |
| `dog_lifespan/` | Dog life tables, the body-size effect, dog-to-human age, trial sizes |
| `damage_repair/` | Models of mechanisms: senescent cells, mutations, clocks, cancer, resilience |
| `interventions/` | Evidence tables, hazard-ratio-to-years math, ranked estimate of years gained |
| `escape_velocity/` | The escape-velocity threshold, periodic repair, population effects |
| `theory/` | `verify_theorems.py`, a computer check of every theorem |
| `tests/` | Test suite |

Each part has its code (`*.py`), its numbers (`results.json`), its charts (`figures/`), and the data behind them (`data/`).

## How to run it

Requires Python 3.11 or newer.

```bash
pip install -r requirements.txt
pytest                 # fast checks: all 19 theorems and the headline numbers (about 1 minute)
pytest -m slow         # re-runs every part and confirms the saved results reproduce (about 6 minutes)
./run_all.sh           # regenerate every result and figure
```

To run one part, change into its folder and run its scripts, for example `cd human_mortality && python3 run_all.py`.

## What is data, what is model, and what is judgement

- **Data.** Life tables and trial results transcribed from published sources. Every transcribed file names its source; see `DATA_SOURCES.md`.
- **Fitted.** Mortality parameters estimated from those life tables.
- **Model.** Theorems and simulations. They are exact statements about the models, not about real bodies. Inputs marked `ILLUSTRATIVE` in the code were chosen to show how a mechanism behaves.
- **Judgement.** The ranked treatment estimates in `interventions/` and the bottleneck ranking in `escape_velocity/` run stated assumptions through arithmetic. Inputs marked `[J]` or `JUDGEMENT` in the code are opinions, not measurements.

## Limits

- Everything rests on the observation that adult death risk grows by a fixed percentage each year. That fits ages 30 to 95 in humans and 2 to 17 in dogs. It fails for the young and is disputed above age 105.
- Most of the theorems restate or sharpen known results. The theory paper lists which are classical and which may be new. No claim of originality should be made until an expert has checked.
- Several items could not be verified. They are listed in `KNOWN_GAPS.md`.
- This work has not been peer reviewed.

## How this was made

This project was produced with Claude (Anthropic's AI model) under human direction. AI agents wrote the code, transcribed the data and drafted the papers. A second set of AI agents then re-ran the code, recomputed the key numbers independently and re-opened the sources. Corrections from that check are listed at the end of the report.

## Licence

Code is under the MIT licence (`LICENSE`). Text, figures and original data files are under CC BY 4.0 (`LICENSE-CC-BY-4.0.md`). Data transcribed from published sources remains subject to the terms of those sources; see `DATA_SOURCES.md`.

## Citing

See `CITATION.cff`.
