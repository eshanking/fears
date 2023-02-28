---
title: 'FEArS: a python package for simulating evolution on arbitrary fitness seascapes'
tags:
  - Python
  - evolution
  - fitness landscape
  - fitness seascape
authors:
  - name: Eshan S. King
    orcid: 0000-0002-0345-3780
    affiliation: 1
  - name: Davis T. Weaver
    orcid: 0000-0003-3086-497X
    affiliation: 1
  - name: Jacob G. Scott
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 2
    orcid: 0000-0003-2971-7673
affiliations:
 - name: Systems Biology and Bioinformatics Program, Case Western Reserve University School of Medicine, USA
   index: 1
 - name: Translational Hematology Oncology Research, Cleveland Clinic Lerner Research Institute, USA
   index: 2
date: 03 November 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The evolution of drug resistance across kingdoms, including in cancer and 
infectious disease, is governed by the same fundamental laws. Modeling 
evolution with genotype-specific dose response curves, collectively forming a
'fitness seascape', enables simulations that include realistic pharmacokinetic 
constraints, more closely resembling the environmental conditions within a 
patient. FEArS (Fast Evolution on Arbitrary Seascapes) is a python package
that enables simulating evolution with fitness seascapes. FEArS can simulate a 
wide variety of experimental conditions with many arbitrary biological 
parameters. FEArS remains computationally efficient despite being an 
agent-based model, even for very large population sizes. FEArS also contains 
powerful and flexible utilities for data analysis, plotting, and experimental
fitness seascape estimation. 

The two core classes for simulating populations and running experiments are 
Population and Experiment, respectively. The Population class includes all
biological parameters relevant for the evolving population under study in 
addition to methods for simulating evolution. The Experiment class includes
parameters for running an experiment with multiple populations, including 
varying pharmacokinetic parameters, number of simulations, and results saving
methods.

## A hybrid agent-based algorithm

FEArS achieves fast runtimes while simulating large populations of evolving
agents by employing what we term a 'hybrid agent-based' approach. When 
possible, populations are stored as vectors of cell numbers $\hat n$, where each 
position in the vector corresponds to a genotype and the number at that 
position gives the number of cells of that type in the population. Then, 
stochastic events such as cell division and cell death are simply drawn from 
poission distributions:

\begin{equation}\label{eq:cell_death}
  \hat n_{d} \sim poisson(r_{d}*\hat n),
\end{equation}

where $\hat n_{d}$ refers to the vector of dead cells of each type, for example.
However, in the mutation step, FEArS switches to a strictly agent-based process.
Here, every mutating agent is enumerated in a vector, where each entry in the 
vector represents the genotype of the agent. Then, mutating agents are randomly
allocated to adjacent genotypes (Fig. \autoref{fig:flowchart}). Since the number
of mutating cells is much smaller than the total population size (i.e., with a 
mutation rate on the order of $10^{-6}$ to $10^{-9}$ per base pair), this 
agent-based step does not compromise computational efficiency.

![FEArS algorithm flow chart. The blue dashed box indicates the portion of the algorithm that is strictly agent-based.\label{fig:flowchart}](fears_flow_chart.png){ width=70% }

By modeling realistic population sizes, FEArS enables investigation of 
population extinction and stochastic evolutionary rescue.

## A suite of useful utilies

In addition to the core population and experiment classes, FEArS includes a 
wide suite of utilities to assist with computational experiments, empirical
data analysis, and results visualization.

- plotter: a broad and flexible plotting utility, including functions for plotting
evolutionary dynamics timetraces, Kaplan-Meier curves, and fitness landscapes.

- pharm: functions for pharmacokinetics, including generating arbitrary 
pharmacokinetic curves and simulating patient nonadherence.

- fitness: functions for pharmacodynamics, including computing fitness landscapes
and fitness seascapes.

- AutoRate: classes and methods for estimating fitness seascapes from 
experimental data

# Statement of need

FEArS enables stochastic simulations of clonally evolving systems 
subject to arbitrary drug concentrations over time. By using an agent-based
algorithm, we are able to model mutation and selection, with evolution arising 
as an emergent phenomena. Furthermore, by allowing for arbitrary population 
sizes, FEArS can model population extinction. Arbitrary population sizes allows
us to simulate how a disease population within a patient may respond to 
therapy. In addition, FEArS models genotype-sprecific dose-response curves, 
allowing for more fine-grained prediction of evolution.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Acknowledgements

<!-- We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project. -->

# References