---
output: 
  html_document: 
    toc: yes
---

This folder contains .RDS objects for the main results of this analysis.

You can load these .RDS files from here instead of from the cache to re-render the figures and tables, or to explore further.

Item descriptions:

##### all_sims.RDS

Estimated total biomass and total energy use for every route in every year sampled, under the actual and null model scenarios.

`timeperiod`: Year.

`abundance`: Total abundance of all species on that route on that year.

`total_biomass`: Total biomass of all individuals on that route in that year for this scenario (null or actual).

`total_energy`: Total metabolic rate of all individuals on that route in that year for this scenario (null or actual).

`countrynum`: Country number identifier, inherited from the raw BBS data.

`statenum`: State number identifier, inherited from the raw BBS data.

`route`: Route identifier, inherited from the raw BBS data.

`routename`: Route name, inherited from the raw BBS data.

`active`: Active or no, inherited from the raw BBS data.

`latitude`: Route latitutde, inherited from the raw BBS data.

`longitude`: Route longitude, inherited from the raw BBS data.

`stratrum`: Stratum, inherited from the raw BBS data.

`bcr`: Bird conservation region, inherited from the raw BBS data.

`routetypeid`: inherited from the raw BBS data.

`routetypedetailid`: inherited from the raw BBS data.

`regionname`: inherited from the raw BBS data.

`source`: "real" or "sim". "Real" is the actual dynamics (incorporating size change); "sim" is the null model dynamics.

`fsource`: `source` as a factor; not used.

`matssname`: The name of the dataset in the cache. Of the form "bbs_rtrg_X\_Y" where X is the route number and Y is the region number.

`simtype`: always "actual"; not used.

##### all_aics_b.RDS

Model identifiers and evaluation results for (generalized) linear models fit to long-term trends in *biomass*.

`model_AIC`: Model AIC score for each model.

`model_AICc`: Model AICc score for each model.

`model_family`: Gaussian or Gamma.

`model_link`: identity or log, for Gaussian or Gamma models, respectively.

`model_formula`: model formula. `timeperiod` refers to `timeperiod` from `all_sims`; `source` refers to `source` in `all_sims`.

`matssname`: dataset identifier, links to `all_sims`.

`model_p`: model p-value in an ANOVA comparing this model to the next-most-complex one (ranked by `modelcomplexity`).

`modelcomplexity`: 1-4, ranking of model complexity from least (1, intercept-only) to most (4, full interaction).

##### all_aics_e.RDS

Model identifiers and evaluation results for (generalized) linear models fit to long-term trends in *energy use*.

All column names are the same as for `all_aics_b`.

##### all_preds_b.RDS

Estimated predicted values and ratio of end:begin values for models fit to *biomass*.

`begin_real`: Fitted value for total biomass at the first time step, taking into account size change.

`end_real`: Fitted value for total biomass at the last time step, taking into account size change.

`begin_sim`: Fitted value for total biomass at the first time step, *without* taking into account size change.

`end_sim`: Fitted value for total biomass at the last time step, *without* taking into account size change.

`ratio_sim`: Ratio of `end_sim` to `begin_sim`; values \>1 indicate a positive slope.

`ratio_real`: Ratio of `end_real` to `begin_real`; values \>1 indicate a positive slope.

`model_cor`: Correlation coefficient of predicted:actual model values; not used.

`model_cor_squared`: `model_cor_squared`; not used.

`model_family`: Gaussian or Gamma; links to `all_aics_b`.

`model_link`: Identity or log; links to `all_aics_b`.

`model_formula`: Links to `all_aics_b`.

`model_begin`: Year of first time step in model.

`model_end`: Year of last time step in model.

`model_span`: `model_end` - `model_begin`; the number of years spanned by the timeseries (Usually 30).

`matssname`: Dataset identifier; links to `all_aics_b` and `all_sims`.

`model_observations`: Number of years sampled in the model. Not necessarily the same as `model_span`, if some years were not sampled.

##### all_preds_e

Estimated predicted values and ratio of end:begin values for models fit to *energy use*.

All columns are the same as for `all_preds_b`, except for energy use rather than biomass.

##### all_cs_compares

Comparisons of community structure (ISD overlap, change in mean body size, and species turnover) from the first 5 to last 5 years of sampling for each route.

`real_overlap`: Overlap score (0-1, with 1 being high overlap) for the ISD for the first 5 compared to last 5 years, *taking into account size change*.

`sim_overlap`: Overlap score, *not taking into account size change*.

`real_turnover`: 1 - `real_overlap`.

`sim_turnover`: 1 - `sim_turnover`.

`begin_years`: Character string naming the years included in the "begin" time chunk.

`end_years`: Character string naming the years included in the "end" time chunk.

`real_begin_mean_mass`: Mean body size of all individuals (pooled) observed in the first 5 year period.

`real_end_mean_mass`: Mean body size of all individuals (pooled) observed in the last 5 year period.

`sim_begin_mean_mass`: Mean body size of all individuals (pooled) observed in the first 5 year period for the *null model*.

`sim_end_mean_mass`: Mean body size of all individuals (pooled) observed in the last 5 year period for the *null model*.

`real_mass_ratio`: Ratio of `real_end_mean_mass` to `real_begin_mean_mass`. Values \>1 indicate an increase in community-wide mean body size over time.

`sim_mass_ratio`: Ratio of `sim_end_mean_mass` to `sim_begin_mean_mass`. Close to 1.

`sp_turnover_bcd`: Bray-Curtis dissimilarity comparing the first 5 year period to the last 5 year period.

`countrynum`: Country number identifier, inherited from the raw BBS data.

`statenum`: State number identifier, inherited from the raw BBS data.

`route`: Route identifier, inherited from the raw BBS data.

`routename`: Route name, inherited from the raw BBS data.

`active`: Active or no, inherited from the raw BBS data.

`latitude`: Route latitutde, inherited from the raw BBS data.

`longitude`: Route longitude, inherited from the raw BBS data.

`stratrum`: Stratum, inherited from the raw BBS data.

`bcr`: Bird conservation region, inherited from the raw BBS data.

`routetypeid`: inherited from the raw BBS data.

`routetypedetailid`: inherited from the raw BBS data.

`regionname`: inherited from the raw BBS data.

`matssname`: dataset identifier, links to other data tables.
