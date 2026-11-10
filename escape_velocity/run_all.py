"""Run every script in this stream and merge their outputs into results.json."""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARTS = [("escape_velocity.py", "results_escape_velocity.json"),
         ("repair_model.py", "results_repair.json"),
         ("population_consequences.py", "results_population.json"),
         ("bottleneck_voi.py", "results_voi.json")]
merged = {"_note": "Merged by run_all.py. JUDGEMENT/ILLUSTRATIVE inputs are labelled in each script and in findings.md."}
for script, res in PARTS:
    subprocess.run([sys.executable, os.path.join(HERE, script)], cwd=HERE, check=True,
                   stdout=subprocess.DEVNULL)
    with open(os.path.join(HERE, res)) as fh:
        merged[script.replace(".py", "")] = json.load(fh)
with open(os.path.join(HERE, "results.json"), "w") as fh:
    json.dump(merged, fh, indent=2, default=float)
print("wrote results.json with sections:", list(merged))
