"""Part 2-4: transcribed evidence tables (mouse ITP, primates, human, dog) and the ITP figure.

Every number below was read from the source URL listed next to it in this session (Nov 2026).
Numbers marked 'unverified' could not be checked against an opened source.
Run:  python3 evidence_tables.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import DATA, FIG, update_results

NIA = "https://www.nia.nih.gov/research/dab/interventions-testing-program-itp/supported-interventions"
# ---------------------------------------------------------------- mouse ITP (UM-HET3), % change in MEDIAN lifespan
itp = [
    # compound, start_age_months, male_pct, female_pct, metric, sig_note, source
    ("Rapamycin (42 ppm, highest dose)", 9, 23, 26, "median", "both sexes", "Miller et al. 2014 Aging Cell; abstract read via https://pubmed.ncbi.nlm.nih.gov/24341993/"),
    ("Rapamycin (from 20 months)", 20, 9, 14, "age at 90% mortality (not median)", "both sexes", "Harrison et al. 2009 Nature; abstract via https://scholars.uthscsa.edu/en/publications/rapamycin-fed-late-in-life-extends-lifespan-in-genetically-hetero/"),
    ("Rapamycin + acarbose", 9, 34, 28, "median", "both sexes", "Strong et al. 2022 Aging Cell e13724; via https://www.lifespan.io/news/combination-of-rapamycin-and-acarbose-extends-lifespan/"),
    ("Rapamycin + metformin", 9, 23, 23, "median", "both sexes", "Strong et al. 2016 Aging Cell, Tables 1-2; https://ueaeprints.uea.ac.uk/58877/4/Manuscript.pdf"),
    ("Acarbose", 4, 22, 5, "median", "males (F: max only)", "Harrison et al. 2014 Aging Cell; https://pubmed.ncbi.nlm.nih.gov/24245565/"),
    ("Acarbose (from 16 months)", 16, 6, 2, "median", "M median; both sexes 90th pct", "Strong et al. 2016; https://ueaeprints.uea.ac.uk/58877/4/Manuscript.pdf"),
    ("17-alpha-estradiol (14 ppm)", 10, 19, 1, "median", "males only", "Strong et al. 2016; https://ueaeprints.uea.ac.uk/58877/4/Manuscript.pdf"),
    ("17-alpha-estradiol (from 16 months)", 16, 19, np.nan, "median", "males only", "Harrison et al. 2021 Aging Cell; https://mouseion.jax.org/stfb2021/88"),
    ("17-alpha-estradiol (from 20 months)", 20, 11, np.nan, "median", "males only", "Harrison et al. 2021; https://mouseion.jax.org/stfb2021/88"),
    ("Canagliflozin", 7, 14, 0, "median", "males only (F: no benefit)", "Miller et al. 2020 JCI Insight; https://doaj.org/article/8346a43abd5742d2afbdcea95b0aa4d1"),
    ("Glycine", 9, 5, 5, "lifespan (4-6%, sexes pooled in abstract)", "both sexes", "Miller et al. 2019 Aging Cell; https://mouseion.jax.org/stfb2019/117"),
    ("NDGA (2500 ppm)", 9, 10, np.nan, "median", "males only", "Strong et al. 2016; https://ueaeprints.uea.ac.uk/58877/4/Manuscript.pdf"),
    ("Astaxanthin", 12, 12, 0, "median", "males only; NOT replicated in 2022 cohort", "Harrison et al. 2024 GeroScience; https://mouseion.jax.org/stfb2024/24 ; replication null per " + NIA),
    ("Meclizine", 12, 8, 0, "median", "males only; NOT replicated in 2022 cohort", "Harrison et al. 2024; https://mouseion.jax.org/stfb2024/24 ; " + NIA),
    ("Protandim", 10, 7, 3, "median", "males only", "Strong et al. 2016; https://ueaeprints.uea.ac.uk/58877/4/Manuscript.pdf"),
    ("Captopril", 5, np.nan, 5, "median", "females small (4-5%); males indeterminate", "Strong et al. 2022; https://www.lifespan.io/news/combination-of-rapamycin-and-acarbose-extends-lifespan/"),
    # failures
    ("Metformin", 9, 7, 0, "median", "NOT significant (M p=0.35)", "Strong et al. 2016; https://ueaeprints.uea.ac.uk/58877/4/Manuscript.pdf"),
    ("Resveratrol", np.nan, 0, 0, "no increase", "no effect (3 tests)", NIA),
    ("Nicotinamide riboside", np.nan, 0, 0, "no increase", "no effect either sex", "Harrison et al. 2021; https://mouseion.jax.org/stfb2021/88"),
    ("Fisetin (continuous or intermittent)", np.nan, 0, 0, "no increase", "no effect either sex", "Harrison et al. 2024; https://mouseion.jax.org/stfb2024/24"),
    ("Simvastatin, enalapril, curcumin, green tea, fish oil, MitoQ, alpha-ketoglutarate, others", np.nan, 0, 0, "no increase", "no effect", NIA),
]
itp_df = pd.DataFrame(itp, columns=["intervention", "start_age_months", "male_pct", "female_pct", "metric", "note", "source"])
itp_df["start_age_months"] = itp_df["start_age_months"]
# start ages for some entries are approximate (ITP standard start ~4-12 months) and labelled so
approx_start = {"Rapamycin (42 ppm, highest dose)", "Rapamycin + acarbose", "Rapamycin + metformin", "Acarbose",
                "17-alpha-estradiol (14 ppm)", "Glycine", "NDGA (2500 ppm)", "Astaxanthin", "Meclizine", "Protandim",
                "Captopril", "Metformin"}
itp_df["start_age_note"] = ["approximate/unverified start age" if n in approx_start else "from source" for n in itp_df.intervention]
with open(f"{DATA}/itp_mouse_lifespan.csv", "w") as fh:
    fh.write("# NIA Interventions Testing Program, UM-HET3 mice. % change in median lifespan vs controls unless 'metric' says otherwise.\n")
    fh.write(f"# Master list of tested compounds (qualitative): {NIA} (page last updated 13 Apr 2026)\n")
    fh.write("# start_age_months: only rows marked 'from source' were checked; others are approximate (unverified).\n")
    itp_df.to_csv(fh, index=False)

# figure
plot = itp_df[itp_df.metric.str.startswith("median") | itp_df.intervention.str.startswith("Glycine")
              | itp_df.intervention.str.startswith("Rapamycin (from")].copy()
plot = plot.iloc[::-1]
fig, ax = plt.subplots(figsize=(8, 6.5))
yy = np.arange(len(plot))
ax.barh(yy + 0.2, plot.male_pct.fillna(0), height=0.4, color="#2a6f97", label="males")
ax.barh(yy - 0.2, plot.female_pct.fillna(0), height=0.4, color="#e07a5f", label="females")
for i, (_, rw) in enumerate(plot.iterrows()):
    if np.isnan(rw.female_pct):
        ax.text(0.5, i - 0.2, "not tested / not reported", fontsize=6, va="center")
ax.set_yticks(yy); ax.set_yticklabels(plot.intervention, fontsize=7)
ax.set_xlabel("% change in median lifespan vs control (rapamycin-from-20-mo row: age at 90% mortality)")
ax.set_title("NIA Interventions Testing Program, genetically mixed mice")
ax.axvline(0, color="k", lw=.5); ax.legend()
fig.tight_layout(); fig.savefig(f"{FIG}/fig_itp_lifespan.png", dpi=150); plt.close(fig)
plot.to_csv(f"{DATA}/fig_itp_lifespan.csv", index=False)

# ---------------------------------------------------------------- other animal evidence
animal = [
    ("Caloric restriction, rhesus monkey, Univ. Wisconsin (adult onset, ~30%)", "survival improved; CR roughly halved death rate vs ad-lib controls", "https://alzforum.org/node/742831 ; https://ideas.repec.org/a/nat/natcom/v5y2014i1d10.1038_ncomms4557.html", "ESTABLISHED for this cohort"),
    ("Caloric restriction, rhesus monkey, NIA (121 monkeys, 25-30% CR; controls not ad-lib)", "no significant survival gain; young-onset CR no benefit; fewer cancers", "https://www.sciencenews.org/article/low-cal-longevity-questioned ; https://alzforum.org/node/39319", "ESTABLISHED for this cohort"),
    ("Joint UW/NIA re-analysis (2017)", "both: age-related disease ~half the rate; survival benefit with adult/late onset only", "https://alzforum.org/node/742831 ; https://www.thesciencebreaker.org/en/breaks/health-physiology/lower-calorie-intake-allows-monkeys-to-live-long-and-prosper", "CONTESTED (onset age, control diet)"),
    ("Dietary restriction in 960 Diversity Outbred female mice (20%, 40% CR, fasting)", "40% CR: median +9 months vs ad lib; genetics explained more lifespan variance than diet; greater weight loss linked to shorter life", "https://www.news-medical.net/news/20241022/Caloric-restriction-extends-mouse-lifespan-and-health-in-new-genetic-study.aspx", "ESTABLISHED (one large study)"),
    ("Dietary restriction across recombinant inbred mouse strains (Liao 2010)", "lifespan response ranges from extension to shortening by strain (numbers not verified)", "https://vivo.colorado.edu/display/pubid_37673 (metadata only; numbers unverified)", "CONTESTED"),
    ("GH-deficient / GH-receptor knockout mice", "GH knockout reported ~60-70% longer-lived (secondary source); Ames/Snell/GHR-KO specific % unverified", "https://www.fightaging.org/archives/2018/05/a-review-of-growth-hormone-in-aging/?nc", "ESTABLISHED direction; magnitudes unverified"),
    ("Human GH-receptor deficiency (Laron; Ecuador & Brazil cohorts)", "~100 people in Ecuador: no diabetes, 1 cancer; Brazil 2019: no increase in overall lifespan", "https://en.wikipedia.org/wiki/Laron_syndrome", "ESTABLISHED (small cohorts)"),
    ("Telomerase AAV gene therapy in mice (one lab)", "median +24% (start 1 y), +13% (start 2 y); no excess cancer reported", "https://doaj.org/article/9b339a05152245db843dd446a714c398", "single lab, not replicated"),
    ("Klotho overexpression (mice, 2005); single klotho injection in 18 old monkeys (2023)", "mouse: extends lifespan (magnitude unverified); monkey working memory ~45%->~60% correct for >=2 weeks", "https://research.joslin.org/kahnlaboratory/publications/suppression-aging-mice-hormone-klotho ; https://www.freethink.com/science/klotho-working-memory", "C/D"),
    ("Taurine supplementation (one lab, Singh 2023 Science)", "mouse median +10-12%; monkey health markers only; 2025 study: blood taurine does not fall with age in humans", "https://genengnews.com/news/lifespan-in-mice-and-monkeys-improved-by-reversing-taurine-deficiency ; https://www.news-medical.net/news/20250605/Study-challenges-taurine-as-a-biomarker-of-aging.aspx", "CONTESTED; ITP test started 2025"),
]
pd.DataFrame(animal, columns=["intervention", "result", "source", "tag"]).to_csv(f"{DATA}/other_animal_evidence.csv", index=False)

# ---------------------------------------------------------------- human evidence (status as looked up Nov 2026)
human = [
    ("Not smoking", "Persistent smokers lose ~10 y; quitting at 30/40/50/60 gains ~10/9/6/3 y", "cohort (50-y British doctors)", "A* (observational, causal consensus)", "established", "https://news.cancerresearchuk.org/2004/06/22/50-year-tobacco-study-shows-smokers-lose-a-decade-of-life/"),
    ("Exercise", "15 min/day vs inactive: all-cause mortality -14%, ~+3 y (observational, 400k Taiwanese); Generation 100 RCT (n=1567, age ~72, 5 y): no significant mortality difference, controls also active", "cohort + RCT", "B (RCT null on mortality; strong observational)", "established habit; RCT mortality benefit unproven", "https://www.eurekalert.org/news-releases/825516 ; https://sciencenorway.no/elderly-exercise-physical-activity/overly-enthusiastic-70-year-olds-messed-up-study-on-physical-activity/1775873"),
    ("Blood-pressure control (SPRINT, intensive vs standard)", "All-cause death HR 0.73 (0.60-0.90), n=9361 high-risk adults >=50, 3.3 y; more hypotension (2.4 vs 1.4%) and acute kidney injury (4.1 vs 2.5%)", "RCT", "A", "standard care", "https://www.wikijournalclub.org/wiki/SPRINT"),
    ("LDL lowering (statins, CTT 2010)", "All-cause mortality rate ratio 0.90 (0.87-0.93) per 1 mmol/L LDL reduction; 26 RCTs, about 170,000 people (the 2005 analysis of 14 RCTs gave 0.88, 0.84-0.91)", "RCT meta-analysis", "A", "standard care", "https://research-portal.uea.ac.uk/en/publications/efficacy-and-safety-of-more-intensive-lowering-of-ldl-cholesterol/"),
    ("Shingles vaccination", "~20% lower 7-year dementia risk (natural experiment, Wales, >280,000 people)", "quasi-experimental", "B", "routine vaccine", "https://www.sciencedaily.com/releases/2025/04/250402122149.htm"),
    ("Caloric restriction (CALERIE)", "220 non-obese adults, 2 y; 25% prescribed, 11.9% achieved; DunedinPACE slowed 2-3% (d=-0.25 at 24 mo); PhenoAge/GrimAge unchanged; authors extrapolate 10-15% lower mortality risk (not observed)", "RCT, biomarker", "B", "no outcome trial", "https://bpb-us-e1.wpmucdn.com/sites.psu.edu/dist/0/9948/files/2023/02/Waziry-Ryan-et-al.-2023-Nature-Aging.pdf"),
    ("Rapamycin, low intermittent (PEARL)", "48 wk, n=114 analysed, ages 50-85: primary endpoint (visceral fat) null; some sex-specific secondary gains (lean mass, pain in women); GI events more common; one anaemia case", "RCT, biomarker", "B", "no outcome trial; off-label use only", "https://static-site-aging-prod.impactaging.com/article/206235/text"),
    ("TORC1 inhibitor RTB101 (immune ageing)", "Phase 3 PROTECTOR 1 (n=1024, >=65): no reduction in respiratory illness (OR 1.07, p=0.65) after earlier positive phase 2", "RCT", "A (negative)", "programme stopped 2019", "https://www.biospace.com/restorbio-announces-that-the-phase-3-protector-1-trial-of-rtb101-in-clinically-symptomatic-respiratory-illness-did-not-meet-the-primary-endpoint"),
    ("Metformin / TAME", "ITP mice: no significant lifespan effect; MeMeMe trial: no effect on cancer, CVD or mortality (secondary source); TAME (planned ~3000 people aged 65-79, composite of death/CVD/cancer/cognitive decline) has no ClinicalTrials.gov record and has not enrolled as of Sept 2026", "RCT (non-ageing) + observational", "C (for ageing)", "TAME not started", "https://www.longevitysciencedaily.com/articles/metformin-tame-trial-never-started ; https://www.fightaging.org/archives/2026/08/a-review-of-the-present-state-of-data-for-metformin-as-a-geroprotective-drug/"),
    ("SGLT2 inhibitor (EMPA-REG OUTCOME)", "Type 2 diabetes + CVD, n=7020, 3.1 y: all-cause death HR 0.68 (0.57-0.82); genital infections 6.4 vs 1.8%", "RCT", "A (in diabetes/CVD) / B (ageing)", "approved for diabetes, heart failure, kidney disease", "https://www.wikijournalclub.org/wiki/EMPA-REG_OUTCOME"),
    ("GLP-1 receptor agonist (SELECT, semaglutide)", "Overweight/obese with CVD, no diabetes, n=17,604, 34 mo: all-cause death HR 0.81 (0.71-0.93); 10 vs 2% stopped for GI events", "RCT", "A (in obesity+CVD) / B (ageing)", "approved for obesity, diabetes, CV risk", "https://www.wikijournalclub.org/wiki/SELECT"),
    ("Senolytics (dasatinib + quercetin)", "Phase 2 RCT, 60 postmenopausal women, 20 wk: bone formation up, resorption not reduced; benefit mainly in high-senescence subgroup; ITP: fisetin no mouse lifespan effect", "RCT, biomarker", "C", "experimental", "https://newsnetwork.mayoclinic.org/discussion/drugs-that-kill-zombie-cells-may-benefit-some-older-women-but-not-all-mayo-clinic-study-finds/ ; https://mouseion.jax.org/stfb2024/24"),
    ("NAD precursors (NR, NMN)", "Meta-analysis of 10 placebo-controlled studies: no improvement in muscle mass or function; ITP: NR no mouse lifespan effect", "RCT meta-analysis", "C (negative)", "supplements; no ageing indication", "https://www.nutraingredients.com/Article/2025/06/25/can-nmn-and-rn-improve-muscle-mass-and-function/"),
    ("Taurine", "Mouse +10-12% (one lab); human blood taurine does not reliably decline with age (2025)", "animal + observational", "C", "ITP test in progress", "https://genengnews.com/news/lifespan-in-mice-and-monkeys-improved-by-reversing-taurine-deficiency"),
    ("Therapeutic plasma exchange", "Aging Cell 2025, 42 people >50: TPE+IVIG -2.61 y and TPE -1.32 y on multi-omic clocks at mid-trial; per a 2026 secondary report, no significant difference vs sham at end of trial (not checked against paper)", "RCT (small), biomarker", "C", "experimental", "https://www.buckinstitute.org/news/clinical-trial-and-multi-omics-analysis-demonstrates-the-impact-of-therapeutic-plasma-exchange-on-biological-age/ ; https://www.pulse.bot/longevity/news/tpe-long-term-effects-in-healthy-elderly-same-as-sham-ff70b699-f048-4462-8670-08889245c03e/"),
    ("Partial epigenetic reprogramming (Life Biosciences ER-100, OSK)", "First participant dosed June 2026 in Phase 1 safety trial in glaucoma and NAION (NCT07290244); no efficacy data", "Phase 1 (eye)", "D (for ageing)", "Phase 1 ongoing", "https://www.biospace.com/press-releases/life-biosciences-announces-first-patient-dosed-in-phase-1-trial-of-er-100-for-optic-neuropathies"),
    ("Gene therapy: telomerase / klotho / follistatin", "Mouse telomerase +13-24% (one lab); klotho: mouse transgenic + monkey cognition; follistatin: no lifespan data found; no controlled human outcome data found", "animal", "C/D", "experimental", "https://doaj.org/article/9b339a05152245db843dd446a714c398"),
]
hdf = pd.DataFrame(human, columns=["intervention", "best_estimate", "evidence_type", "grade", "status_nov2026", "source"])
hdf.to_csv(f"{DATA}/human_evidence.csv", index=False)

# ---------------------------------------------------------------- dog evidence
dog = [
    ("Diet restriction (Purina Labrador lifetime study)", "48 Labradors, paired littermates, 25% less food from puppyhood: median 13.0 vs 11.2 y (+1.8 y, ~15%); later signs of ageing", "randomised paired, lifetime", "A (single small study)", "https://newscenter.purina.com/LifeSpanStudy"),
    ("Being overweight (observational)", "50,787 dogs, 12 breeds: overweight dogs lived 5 months to 2.5 y less depending on breed", "retrospective cohort", "B (observational)", "https://news.liverpool.ac.uk/2019/01/03/research-reveals-overweight-dogs-may-live-shorter-lives/"),
    ("Dental cleaning (observational)", "2.37 M Banfield dogs: each extra yearly cleaning HR 0.817 for death; authors warn not causal", "retrospective cohort", "B (observational, confounding likely)", "https://jaaha.kglmeridian.com/view/journals/aaha/55/3/article-p130.xml"),
    ("Neutering (observational)", "Banfield: median ~30 weeks longer spayed vs intact females, ~3 weeks for males; >40,000 dogs: longer life but more cancer deaths", "retrospective cohort", "B (observational, CONTESTED)", "https://jaaha.kglmeridian.com/view/journals/aaha/55/3/article-p130.xml ; https://doaj.org/article/6baa7114fc9e4b0592ea0bc333220527"),
    ("Rapamycin, small trials", "10-wk trial: better heart-function indices; 6-month RCT of 17 dogs (2023): no significant echo differences; mild GI side effects", "small RCTs", "B-/C", "https://www.frontiersin.org/journals/veterinary-science/articles/10.3389/fvets.2023.1168711/full"),
    ("Rapamycin, TRIAD (Dog Aging Project)", "Planned 580 dogs >=7 y, 20-55 kg, 1:1 weekly rapamycin vs placebo for 12 mo + 24 mo follow-up; primary endpoint survival; >180 enrolled Feb 2025; planned end Nov 2029; no results", "RCT ongoing", "pending", "https://www.avma.org/news/7m-grant-rescues-dog-study-investigating-rapamycin-canine-aging ; https://jadams5.github.io/aging_wiki/studies/coleman-2025-triad-protocol"),
    ("LOY-001 (Loyal)", "Injectable that lowers IGF-1 in large/giant dogs; FDA accepted RXE (effectiveness) section Nov 2023; safety and manufacturing sections not reported as complete; not approved", "company data (not public)", "C/D", "https://www.dvm360.com/view/fda-determines-drug-for-lifespan-extension-in-large-dogs-to-have-a-reasonable-expectation-of-effectiveness ; https://loyal.com/downloads/Loyal-product-brief-for-vets-v8.pdf"),
    ("LOY-002 (Loyal)", "Daily pill for dogs >=10 y, >=14 lb, targets age-related metabolic dysfunction (ingredients confidential); RXE accepted Feb 2025, safety (TAS) accepted (announced 13 Jan 2026); manufacturing (CMC) outstanding in latest sources; conditional approval not granted as of last dated source found (22 Sep 2026); STAY trial 1,300 dogs, randomised, placebo-controlled, >=4 y, no lifespan results", "company data + ongoing RCT", "C (pending STAY)", "https://loyal.com/posts/loy-002-tas-announcement ; https://www.glp1pets.com/loy-002 ; https://fortune.com/2026/04/11/dog-lifespan-longevity-loyal-drug-fda/ ; https://longevity.technology/news/loyals-big-dog-lifespan-drug-clears-third-fda-hurdle/"),
    ("LOY-003 (Loyal)", "Daily pill lowering GH/IGF-1 in large/giant dogs (with Crinetics); FDA accepted RXE, announced 16 Sep 2026; not approved", "company data", "C/D", "https://longevity.technology/news/loyals-big-dog-lifespan-drug-clears-third-fda-hurdle/"),
]
ddf = pd.DataFrame(dog, columns=["intervention", "best_estimate", "evidence_type", "grade", "source"])
ddf.to_csv(f"{DATA}/dog_evidence.csv", index=False)

update_results("evidence_tables", {
    "itp_rows": len(itp_df), "human_rows": len(hdf), "dog_rows": len(ddf),
    "purina_median_gain_years": 13.0 - 11.2,
    "purina_implied_HR_first_order_b0.25": float(np.exp(-0.25 * 1.8)),
    "rapa_acarbose_combo_vs_singles_male_pct": {"combo": 34, "rapa_42ppm": 23, "acarbose": 22,
                                                 "note": "different cohorts and doses"},
})
print(itp_df[["intervention", "male_pct", "female_pct"]].to_string())
