# Robust inference methods of diagnostic test accuracy meta-analysis for influential outlying studies via density power divergence.

## Main script

The central file is:

- `src/application_example.R`

Run this script to execute the full analysis workflow.

The proposed robust DTAMA method is implemented in:

- `src/functions/robustDTAMA.R`

## Included files

- `data/MMSE.csv`
  - The dataset from: Arevalo‐Rodriguez I, Smailagic N, Roqué i Figuls M, et al. Mini‐Mental State Examination (MMSE) for the detection of Alzheimer’s disease and other dementias in people with mild cognitive impairment (MCI). Cochrane Database Syst Rev. 2015;2015(3):CD010783. doi:10.1002/14651858.CD010783.pub2
- `application_example.R`
- `src/functions/`

## Requirements

- R (implemented and tested with R 4.5.1)
- R packages:
  - `faraway`
  - `mada`
  - `patchwork`
  - `ggplot2`

If needed, install packages manually in R:

```r
install.packages(c("faraway", "mada", "patchwork", "ggplot2"))
```

## How to run

Set your working directory to the project folder, then run:

```r
source("src/application_example.R")
```

## Output

The script writes results and figures to:

- `output/application_example/`

When successful, it prints:

`All analyses completed!`
