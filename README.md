

This repo contains code to reproduce the analyses in "Temporal changes in the individual size distribution modulate the long-term trends of biomass, and energy use of North American breeding bird communities", by Renata Diaz and S. K. Morgan Ernest, in press at _Global Ecology and Biogeography_. 

This code depends on the following R packages. If you do not have them installed, you can install them by running `install.packages("<PACKAGE_NAME>")` in R. 

```
dplyr
tidyr
vegan
here
maps
mclust
MuMIn
R.utils
```

To re-run all analyses:

- Download this repo (You can download it as a .zip file by clicking on the green "Code" button and selecting "Download as .zip").
- Open the included .Rproj file.
- Run the script `analysis/pipelines/pipeline.R`. Results files will be stored in the `results` directory. Note that, by default, this script will only run analyses for the first 200 (of 1331 possible) sites. This takes about half an hour on a MacBook Air with 16GB of memory. Running for all sites typically takes about 3 hours on the same machine, and must be allowed to run in the same R session for the entire analysis.

You can then load the individual .RDS files (which are described in detail in `analysis/results/object_descriptions.md` to explore them.

Or, you can reproduce the figures and tables included in the manuscript and appendices by opening and rendering the RMarkdown documents in `analysis/writing/final-version/source_documents`.
