# Data folder

This folder is intentionally empty in the public repository. The Lifelines
Cohort Study consent does not allow individual-level data to be shared in a
public archive. Bona-fide researchers can apply for access via the Lifelines
Research Office (<https://www.lifelines.nl/researcher>), after which the
dataset can be opened inside the Lifelines secure HPC environment.

## Files expected by the analysis scripts

The scripts expect the following files at the paths configured in `01_Config/01_setup.R`:

| File | Used by | Description |
|---|---|---|
| `Relatives_with_project_id_M_and_F` | 02 | Pedigree (PSEUDOIDEXT, FatherID, MotherID). |
| `HouseID_Final_with_project_id`     | 02, 04 – 05 | Household identifier — models the shared environment (*c*). |
| `data_for_baseline.Rdata`           | 02, 03 | Baseline phenotype table with CTQ-SF, Loneliness, SPF-IL, LTE / LDI sum scores (the output of `02a_score_construction.R` or an equivalent earlier extraction). |
| `data_for_baseline_transformed.Rdata` | 05 | Same table after rank-based inverse normal transformation (Suppl. Table S10). |
| `adult_personality_OUTPUT.csv`      | 02 | Lifelines NEO release. The three analytic columns are: Neuroticism = NEO-PI-R 48-item six-facet domain sum (range 48 – 240); Extraversion = NEO-FFI domain SUM3 (11 items administered, rescaled to a 12 – 60 analytic scale via the Lifelines `(SUM2/NVALID) × 12` rule); Conscientiousness = NEO-FFI domain SUM3 (12 items, range 12 – 60). All three personality variables are z-standardised before variance decomposition. See the Lifelines NEO wiki for full scoring details. |
| `phenotypes_lifelines.Rdata`        | 02 (optional) | Output of `02a_score_construction.R`. Skip 02a if this file is already present. |

## Raw item-level files (used only by `02a_score_construction.R`)

If you want to re-run the item-level cleaning from scratch, the following raw
Lifelines distributions are needed under `data/raw/`:

| File | Contents |
|---|---|
| `CTQ_items.csv` | 28 CTQ-SF items per participant (`ctq_*_adu_q_01` … `_28`). |
| `Loneliness_items.csv` | 6 De Jong Gierveld items (`loneliness_*_adu_q_01` … `_06`). |
| `SocialSupport_items.csv` | 9 SPF-IL items (`spfil_socsup_adu_q_01` … `_09`). |
| `Mental_Health_items.csv` | 12 BRUGHA + 12 GLLM items per assessment wave. |

## Public item-level documentation (Lifelines Wiki)

Every item code is documented (with original wording, response options, and
per-wave administration schedule) on the **Lifelines Wiki**:

| Instrument | Wiki page |
|---|---|
| Childhood Trauma Questionnaire (CTQ-SF) | <https://wiki.lifelines.nl/doku.php?id=childhood_trauma_ctq> |
| De Jong Gierveld Loneliness Scale | <https://wiki.lifelines.nl/doku.php?id=loneliness_ls> |
| SPF-IL Social Production Function Instrument | <https://wiki.lifelines.nl/doku.php?id=social_functioning_spf-il> |
| List of Threatening Experiences (LTE, Brugha) | <https://wiki.lifelines.nl/doku.php?id=stress_lte> |
| Long-term Difficulties Inventory (LDI) | <https://wiki.lifelines.nl/doku.php?id=stress_ldi> |
| NEO Personality Inventory (overview) | <https://wiki.lifelines.nl/doku.php?id=personality_neo> |
| NEO sum-score documentation (SUM1 / SUM2 / SUM3) | <https://wiki.lifelines.nl/doku.php?id=personality_sum_scores_neo> |
