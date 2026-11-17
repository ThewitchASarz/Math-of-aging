"""Experiment 014, second part: more trials for the cases that sat on the edge of a prediction band.

All with a treatment that does nothing, 2,500 simulated trials each, so that the slope test's false-alarm
rate is known to about half a percentage point. Also one diagnostic: sick dogs withdrawn at the moment of
death (window 0) instead of up to 6 months before, to see whether the window matters.

Run:  python3 borderline.py            (about 5 minutes on 2 cores)
      python3 borderline.py --quick
"""
import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from multiprocessing import Pool

import dropout as d

SEED = 20261119
EXTRA = {
    "reference": dict(),
    "sick dogs leave, equally": dict(sick=(0.30, 0.30)),
    "sick dogs leave, equally, at the moment of death": dict(sick=(0.30, 0.30), window=0.0),   # diagnostic
    "frail treated dogs leave": dict(frail=0.05),
}


def main():
    quick = "--quick" in sys.argv
    nsim = 20 if quick else 2500
    d.PATTERNS.update(EXTRA)
    a, b_class, shares, _ = d.br.size_classes()
    v = d.rt.frailty_variance()
    cases = [(p, "no effect", nsim, a, b_class, shares, v, SEED + 1000 * i) for i, p in enumerate(EXTRA)]
    with Pool(2, initializer=d.PATTERNS.update, initargs=(EXTRA,)) as pool:
        rows = pool.map(d.one_case, cases, chunksize=1)
    out = dict(note="More trials for the edge cases of experiment 014. Treatment does nothing in every case. "
                    "Same trial and analysis as dropout.py. All dropout patterns are ASSUMPTIONS.",
               seed=SEED, trials_per_case=nsim, cases=rows)
    with open(os.path.join(d.HERE, "borderline_quick.json" if quick else "borderline.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in rows:
        print(f"{r['pattern']:50s} | slope {r['slope_test_fires']:.4f} "
              f"(slow {r['slope_test_fires_pointing_to_slowed_aging']:.4f}) s-hat {r['median_slope_ratio']:.3f} | "
              f"risk {r['risk_test_fires']:.4f} r-hat {r['median_risk_ratio']:.3f}")


if __name__ == "__main__":
    main()
