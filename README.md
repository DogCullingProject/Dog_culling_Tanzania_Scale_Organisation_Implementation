# Dog culling in Tanzania

Analysis code for:

> "Any dog that was roaming around… we deal with it": investigating the scale, organisation, implementation and motivations behind dog culling in a rabies-endemic setting.

## Contents

| File | What it does |
|---|---|
| `dog_culling_analysis.R` | All quantitative analyses: survey descriptives, logistic regression, power analysis, Figure 1, the national estimate (cull size model, simulation, sensitivity analyses), Figure S1 and Tables 2, 3 and S2–S5 |
| `code/functions/trunc_dist_funs.R` | Zero-truncated negative binomial and Poisson functions used by the main script |

The script is split into numbered sections. The header lists which section produces each table and figure.

## Data

The survey data are not included, as they could identify the local government authorities (LGAs) and officers who took part. They are available on request from the corresponding author.

The script expects three files in `data_shared/`:

- `occurrence.csv`: survey of occurrence (49 LGAs)
- `practice.csv`: survey of practice (16 LGAs)
- `councils_2022.csv`: all LGAs in the 2022 census (194 units)

Table S1, the funding sources in Table 3 and the number of regions with culling were produced from data that include LGA names, so they are not reproduced here.

## Running the code

Open R in the project folder and run `dog_culling_analysis.R` from the top. Section 4 needs a working Stan installation for brms (see https://mc-stan.org/cmdstanr/ or the rstan installation guide). Section 4d takes about 15 minutes; its results are saved to `output/table_s5/` as they finish, so it can be stopped and restarted.

R packages: brms, readr, dplyr, stringr, tibble, fitdistrplus, truncdist, logistf, ggplot2, patchwork, DHARMa, MASS, survival, posterior. Package versions are printed at the end of the script (`sessionInfo()`).
