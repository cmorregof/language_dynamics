# Grammar Competition and Contact-Induced Instability in the Old-to-Middle English Transition

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains reproducible simulation code and generated figures for:

> **Grammar Competition and Contact-Induced Instability in the Old-to-Middle English Transition: A Computational Case Study Using the Symmetric Language Dynamical Equation**
>
> Carlos Manuel Orrego Franco; Juan Carlos Riaño-Rojas
>
> *Computational case study, 2026*

We implement Mitchener's (2003) fully symmetric language dynamical equation to explore how Scandinavian contact could have transiently destabilized competition among Old English grammatical varieties. The model has three historical phases: pre-contact competition among four Old English dialect components, a contact phase that introduces an Old Norse component, and a post-contact phase that returns to four components. The analysis is exploratory and methodological; it is not a causal claim about the Old-to-Middle English transition.

- **Preprint**: to be added upon upload
- **Paper DOI**: to be added upon publication

![No-contact vs moderate contact](figures/lde_comparison_nocontact_vs_moderate.png)

## Scope

This repository provides:

- the main simulation script;
- seven generated PNG figures;
- the LaTeX manuscript source and compiled PDF;
- named exploratory parameter regimes;
- transient displacement diagnostics `D_H` and `D_M2` printed by the script.

It does **not** establish a historical causal mechanism. The simulations are intended as an exploratory computational case study.

## Repository Structure

```text
language_dynamics/
├── LICENSE
├── README.md
├── requirements.txt
├── simulation.py
├── manuscript/
│   ├── grammar_competition_oe_me_case_study.tex
│   └── grammar_competition_oe_me_case_study.pdf
└── figures/
    ├── lde_comparison_nocontact_vs_moderate.png
    ├── lde_heatmap_m2.png
    ├── lde_moderate_diagnostics_entropy.png
    ├── lde_moderate_diagnostics_m2.png
    ├── lde_moderate_frequencies.png
    ├── lde_scenario_entropy.png
    └── lde_scenario_m2.png
```

## Reproducing the Figures

### Requirements

Python 3 with NumPy, SciPy, and Matplotlib.

```bash
pip install -r requirements.txt
```

### Run the simulation

```bash
python3 simulation.py
```

The script writes all figures to `figures/` and prints numerical diagnostics to stdout. It uses SciPy's `RK45` integrator with `rtol=1e-10` and `atol=1e-12`.

## Manuscript

The manuscript source is available in `manuscript/grammar_competition_oe_me_case_study.tex`, with a compiled PDF at `manuscript/grammar_competition_oe_me_case_study.pdf`.

## Scenarios

`q1 / q2 / q3` denote learning fidelity in the pre-contact, contact, and post-contact phases. `a` is average cross-grammar intelligibility. `ON fraction` is the Old Norse fraction introduced at contact onset.

| Scenario | q1 | q2 | q3 | a | ON fraction |
|---|---:|---:|---:|---:|---:|
| Conservative contact | 0.90 | 0.75 | 0.85 | 0.80 | 0.15 |
| Moderate contact | 0.90 | 0.65 | 0.80 | 0.80 | 0.25 |
| Strong contact | 0.90 | 0.55 | 0.75 | 0.75 | 1/3 |
| Canonical Mitchener a=0.5 | 0.90 | 0.60 | 0.80 | 0.50 | 1/3 |
| No-contact baseline | 0.90 | 0.90 | 0.90 | 0.80 | 0.00 |

Initial Old English weights are 0.35 / 0.35 / 0.20 / 0.10 for West Saxon, Mercian, Northumbrian, and Kentish in all scenarios except strong contact, which uses 0.30 / 0.35 / 0.20 / 0.15.

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
  title  = {Grammar Competition and Contact-Induced Instability in the Old-to-Middle English Transition},
  author = {Orrego Franco, Carlos Manuel and Ria{\~n}o-Rojas, Juan Carlos},
  year   = {2026},
  note   = {Computational case study using the symmetric language dynamical equation}
}
```

## License

Code in this repository is released under the [MIT License](LICENSE).

## Contact

- Carlos Manuel Orrego Franco
- ORCID: [0009-0001-9163-5137](https://orcid.org/0009-0001-9163-5137)
- Juan Carlos Riaño-Rojas
- ORCID: [0000-0002-5719-2854](https://orcid.org/0000-0002-5719-2854)
- Universidad Nacional de Colombia
- [GitHub: @cmorregof](https://github.com/cmorregof)
