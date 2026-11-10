"""Run the whole stream in order (fits first; later scripts read results.json)."""
import subprocess
import sys
import os

here = os.path.dirname(os.path.abspath(__file__))
for s in ["fit_mortality_laws.py", "le_sensitivity.py", "why_gompertz.py", "limits_extremes.py",
          "no_aging.py", "summarize.py"]:
    print(f"== {s}", flush=True)
    subprocess.run([sys.executable, s], cwd=here, check=True, stdout=subprocess.DEVNULL)
print("done")
