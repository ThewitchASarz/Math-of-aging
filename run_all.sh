#!/usr/bin/env bash
# Regenerate every result and figure. Takes about six minutes.
set -euo pipefail
cd "$(dirname "$0")"

(cd human_mortality && python3 run_all.py)
(cd dog_lifespan && python3 run_all.py)
(cd damage_repair && for s in evolution saturating_removal clocks cancer resilience; do python3 "$s.py"; done)
# interventions reads the fitted values written by dog_lifespan, so it runs after it
(cd interventions && for s in translation_math evidence_tables monte_carlo trial_design; do python3 "$s.py"; done)
(cd escape_velocity && python3 run_all.py)
(cd theory && python3 verify_theorems.py)
echo "All parts finished."
