# Grammar Competition and Contact-Induced Transient Perturbations

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains reproducible code, generated figures, CSV summaries, and manuscript source for:

> **Grammar Competition and Contact-Induced Transient Perturbations in the Old-to-Middle English Transition: A Computational Case Study Using the Symmetric Language Dynamical Equation**
>
> Carlos Manuel Orrego Franco; Juan Carlos Riaño-Rojas
>
> *Computational case study, 2026*

We implement Mitchener's (2003) fully symmetric language dynamical equation to ask whether Scandinavian contact can be represented, cautiously, as a temporary perturbation in Old English grammar competition. The study is exploratory and methodological: the historically motivated scenarios show measurable transient displacement in entropy and dimension-normalized concentration, but they do not show a different long-run attractor or establish a causal account of the Old-to-Middle English transition.

- **Preprint**: to be added upon upload
- **Paper DOI**: to be added upon publication
- **Repository**: <https://github.com/cmorregof/language_dynamics>

![No-contact vs moderate contact](figures/lde_comparison_nocontact_vs_moderate.png)

## Repository Structure

```text
language_dynamics/
├── README.md
├── requirements.txt
├── simulation.py
├── results/
│   ├── scenario_results.csv
│   ├── ablation_results.csv
│   └── redistribution_sensitivity.csv
├── manuscript/
│   ├── grammar_competition_oe_me_case_study.tex
│   └── grammar_competition_oe_me_case_study.pdf
└── figures/
    ├── lde_moderate_frequencies.png
    ├── lde_moderate_diagnostics_entropy.png
    ├── lde_moderate_diagnostics_concentration.png
    ├── lde_moderate_diagnostics_m2_raw.png
    ├── lde_scenario_entropy.png
    ├── lde_scenario_concentration.png
    ├── lde_scenario_m2_raw.png
    ├── lde_comparison_nocontact_vs_moderate.png
    ├── lde_heatmap_final_c.png
    ├── lde_moderate_ablation_entropy_concentration.png
    └── lde_redistribution_sensitivity.png
```

## Reproducing Figures and Tables

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the complete simulation workflow:

```bash
python3 simulation.py
```

The script regenerates all PNG figures in `figures/` and writes table-ready CSV summaries in `results/`. It uses SciPy `RK45` with `rtol=1e-10` and `atol=1e-12`.

## Scenarios

`q1 / q2 / q3` denote learning fidelity before, during, and after contact. `a` is average cross-grammar intelligibility. `ON fraction` is a contact-intensity parameter, not a literal demographic estimate.

| Scenario | q1 | q2 | q3 | a | ON fraction |
|---|---:|---:|---:|---:|---:|
| Conservative contact | 0.90 | 0.75 | 0.85 | 0.80 | 0.15 |
| Moderate contact | 0.90 | 0.65 | 0.80 | 0.80 | 0.25 |
| Strong contact | 0.90 | 0.55 | 0.75 | 0.75 | 1/3 |
| Canonical Mitchener comparison | 0.90 | 0.60 | 0.80 | 0.50 | 1/3 |
| No-contact baseline | 0.90 | 0.90 | 0.90 | 0.80 | 0.00 |

## Outputs

The main diagnostics are normalized entropy, raw `M2`, dimension-normalized concentration `C_n`, transient displacement integrals `D_H` and `D_C`, and a cautious raw `D_M2` diagnostic. `D_C` is the preferred concentration metric when comparing four- and five-component phases.

## Citation

If you use this code, please cite Mitchener's original model:

```bibtex
@article{mitchener2003bifurcation,
  title   = {Bifurcation analysis of the fully symmetric language dynamical equation},
  author  = {Mitchener, W. Garrett},
  journal = {Journal of Mathematical Biology},
  volume  = {46},
  pages   = {265--285},
  year    = {2003},
  doi     = {10.1007/s00285-002-0172-8}
}
```

Please also cite this repository if it supports your work:

```bibtex
@misc{orrego2026grammarcompetition,
  title  = {Grammar Competition and Contact-Induced Transient Perturbations
            in the Old-to-Middle English Transition},
  author = {Orrego Franco, Carlos Manuel and Ria{\~n}o-Rojas, Juan Carlos},
  year   = {2026},
  note   = {Computational case study using the symmetric language dynamical equation}
}
```

## License

Code and manuscript support files in this repository are released under the [MIT License](LICENSE).

## Contact

Carlos Manuel Orrego Franco, Universidad Nacional de Colombia

ORCID: [0009-0001-9163-5137](https://orcid.org/0009-0001-9163-5137)

GitHub: [@cmorregof](https://github.com/cmorregof)

Juan Carlos Riaño-Rojas, Universidad Nacional de Colombia

ORCID: [0000-0002-5719-2854](https://orcid.org/0000-0002-5719-2854)
