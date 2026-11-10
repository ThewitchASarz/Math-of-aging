"""Run every script of the dog_lifespan stream in dependency order."""
import os
import runpy

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
if os.path.exists("results.json"):
    os.remove("results.json")
for s in ("fit_mortality.py", "size_scaling.py", "dog_years.py", "causes.py", "power.py"):
    print(f"===== {s}")
    runpy.run_path(os.path.join(HERE, s), run_name="__main__")
