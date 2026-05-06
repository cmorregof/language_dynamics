# Grammar Competition and Contact-Induced Instability in the Old-to-Middle English Transition

This repository provides reproducible code for a computational case study using Mitchener (2003)'s fully symmetric language dynamical equation. The script explores how adding an Old Norse grammar during a contact phase can transiently change competition among Old English grammatical varieties. The scenarios are exploratory and methodological: they illustrate model behavior under named parameter regimes and should not be read as a causal claim about the historical transition from Old English to Middle English.

## Reproducing the figures

```bash
python simulation.py
```

Running the script writes seven PNG figures to `figures/`.

## Parameters

| Scenario | q1 | q2 | q3 | a | ON fraction |
|---|---:|---:|---:|---:|---:|
| Canonical Mitchener a=0.5 | 0.90 | 0.60 | 0.80 | 0.50 | 1/3 |
| Conservative contact | 0.90 | 0.75 | 0.85 | 0.80 | 0.15 |
| Moderate contact | 0.90 | 0.65 | 0.80 | 0.80 | 0.25 |
| Strong contact | 0.90 | 0.55 | 0.75 | 0.75 | 1/3 |
| No-contact baseline | 0.90 | 0.90 | 0.90 | 0.80 | 0.00 |

## Citation

Mitchener, W. G. (2003). DOI: [10.1007/s00285-002-0172-8](https://doi.org/10.1007/s00285-002-0172-8).

Author: Carlos Manuel Orrego Franco. License: MIT.
