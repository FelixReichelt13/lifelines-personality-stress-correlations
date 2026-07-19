# 04_Main_Analysis

Primary pedigree-based ACE models (ASReml-R).

| Order | Script | Output |
|-------|--------|--------|
| 1 | [`04_univariate_heritability.R`](04_univariate_heritability.R) | Figure 1; S8 |
| 2 | [`04b_bivariate_stress_stress.R`](04b_bivariate_stress_stress.R) | Figure 2; S13, S15 |
| 3 | [`04c_bivariate_personality_stress.R`](04c_bivariate_personality_stress.R) | Figure 3; S14, S16 |
| 4 | [`04d_significance_audit.R`](04d_significance_audit.R) | Independent audit (`table_audit_*`, `audit_summary.txt`) |

Run scripts 04 and 04b before script 05 (Sensitivity). Script 04d only needs the INTERNAL TSVs from 04–04c.
