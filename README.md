# ms_niche_interests

Code repository for the analysis presented in the manuscript: "Same structure, different stability: The lasting legacy of modelling assumptions".

This project evaluates whether food-web reconstruction frameworks produce interchangeable ecological starting points, or whether the assumptions embedded in their network-generating rules leave a lasting imprint on community dynamics. The analysis integrates food-web generation, bioenergetic simulation, topology analysis, and stability comparison within a common computational workflow.

## Overview

The study compares seven food-web reconstruction approaches spanning structural, statistical, and trait-based philosophies. Each is used to generate an initial interaction network, which is then embedded in a shared bioenergetic dynamical model. The key question is whether initially distinct networks converge structurally under ecological dynamics, and whether any convergence is accompanied by convergence in equilibrium biomass and stability properties.

The repository therefore contains the computational pipeline used to:

- generate candidate food webs under multiple reconstruction frameworks;
- simulate community dynamics under a common bioenergetic formulation;
- quantify post-dynamics network structure and persistence;
- compare local stability and network-level behaviour across models;
- test whether model-specific assumptions explain variation beyond species richness.

## Repository structure

```text
.
├── README.md                     # project overview and analysis workflow
├── index.qmd                     # manuscript source (Quarto)
├── references.bib                # reference database
├── _quarto.yml                  # Quarto project configuration
├── requirements.txt             # Python environment for document tooling
├── code/
│   ├── Project.toml             # Julia project and dependency manifest
│   ├── code.Rproj               # R project file
│   ├── 01_BuildNetworks.jl      # generate synthetic food webs under each model
│   ├── 02_DynamicSimulation.jl  # simulate bioenergetic dynamics and compute stability metrics
│   ├── 03_EcologicalRealism.jl  # assess network realism and filtering criteria
│   ├── 04_Topology.jl           # compute topology statistics before/after dynamics
│   ├── 05_EcologicalRealismAnalysis.R
│   ├── 06_TopologyAnalysis.R    # phenotypic trajectory analysis and convergence tests
│   ├── 08_StabilityAnalysis.R   # model comparison and stability summaries
│   ├── 09_Robustness.jl         # robustness analyses
│   ├── 10_MultivariateDistance.R
│   ├── 11_SEM.R                 # structural equation modelling
│   ├── lib/                     # shared Julia/R utility functions
│   ├── networks/                # generated network objects and saved data
│   └── outputs/                 # CSV outputs used for figures and analyses
├── tables/                       # tabular results used in the manuscript
├── figures/                     # publication figures
├── docs/                        # rendered manuscript outputs
├── supplement/                  # supplementary material
├── LICENSE
└── .github/
```

## Scientific workflow

### 1. Food-web generation

`code/01_BuildNetworks.jl` constructs synthetic food webs under seven generative models, including structural models (Niche, Cascade, Random, MaxEnt) and trait-informed models (Latent Trait Model, ATN, ADBM). Each replicate fixes a common species pool and parameter restrictions so that differences among networks can be attributed to reconstruction assumptions rather than arbitrary variation in richness or dynamical setup.

The script stores adjacency matrices, body masses, and metabolic classes in the generated network objects.

### 2. Dynamic simulation

`code/02_DynamicSimulation.jl` embeds each initial network in a shared bioenergetic food-web model implemented with `EcologicalNetworksDynamics.jl`. The pipeline includes:

- burn-in simulation to identify non-persistent species or biologically implausible links;
- filtering of extinct or unrealistic components;
- re-simulation to a steady state for the realised network;
- calculation of persistence, biomass diversity, flux inequality, interaction skewness, resilience, and reactivity.

This stage is central to the study, because it separates the effect of network assembly from the effect of ecological dynamics under identical dynamical rules.

### 3. Topological and realism assessment

The topology and realism analyses examine whether the ecological dynamics reshape network structure in a model-independent way.

- `code/03_EcologicalRealism.jl` evaluates structural realism, including disconnected and illogical configurations.
- `code/04_Topology.jl` extracts network metrics from the initial and realised food webs.
- `code/05_EcologicalRealismAnalysis.R` visualises changes in realism before and after simulation.
- `code/06_TopologyAnalysis.R` performs phenotypic trajectory analysis (PTA), comparing trajectories in a multivariate topology space and quantifying convergence toward a common realised topology.

These analyses test whether food-web models start in distinct parts of topology space but collapse toward a shared realised region through extinction-driven pruning.

### 4. Stability comparisons

`code/08_StabilityAnalysis.R` is the main comparative analysis of dynamical outcomes. It evaluates differences among reconstruction frameworks for multiple responses, including persistence, biomass diversity, and local stability metrics. The workflow includes effect-size decomposition and adjusted comparisons that account for equilibrium species richness, enabling the study to distinguish model effects from richness effects.

### 5. Multivariate and causal analysis

The remaining scripts relate network structure to ecological behaviour beyond simple pairwise comparisons:

- `10_MultivariateDistance.R` links topological distance to dynamical distance in multivariate space;
- `11_SEM.R` fits multi-group structural equation models to separate direct and indirect pathways linking richness, topology, and stability-related intermediates.

Together, these analyses identify which aspects of ecosystem behaviour are shaped by reconstruction assumptions and which are explained primarily by realised richness.

## Key analyses and outputs

The analysis produces a suite of result tables and figures, including:

- topology trajectories and centroid convergence summaries;
- model comparisons of persistence, resilience, and reactivity;
- decomposition of variance by model and realised richness;
- multivariate distance summaries linking topology and stability;
- SEM parameter estimates and model-path comparisons.

Outputs are stored in `code/outputs/`, `tables/`, and `figures/`, and are used directly in the manuscript and supplementary material.

## Reproducibility

### Julia

From the repository root:

```bash
cd code
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then run the primary simulation steps:

```bash
julia --project=. 01_BuildNetworks.jl
julia --project=. 02_DynamicSimulation.jl
julia --project=. 09_Robustness.jl
```

The network-generation and simulation scripts are the core computational steps for reproducing the study data.

### R analyses

From the `code/` directory:

```bash
Rscript 05_EcologicalRealismAnalysis.R
Rscript 06_TopologyAnalysis.R
Rscript 08_StabilityAnalysis.R
Rscript 10_MultivariateDistance.R
Rscript 11_SEM.R
```

These scripts reproduce the statistical and visualization components of the manuscript.

## Interpretation

The central result of the project is that structural convergence does not necessarily imply dynamical equivalence. Different food-web reconstruction frameworks may generate distinct starting networks, yet ecological dynamics can prune these networks toward a common realised structure. However, the equilibrium biomass and stability properties remain model-dependent, indicating that assumptions about how interactions are generated continue to shape ecological inference even after filtering by persistence.

This repository therefore supports a publication-level investigation of how reconstruction assumptions affect ecological conclusions about food-web stability and structure.

## License

See the project `LICENSE` file for licensing information.
