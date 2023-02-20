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

# Statement of need

FEArS enables stochastic simulations of clonally evolving systems 
subject to arbitrary drug concentrations over time. By using an agent-based
algorithm, we are able to model mutation and selection, with evolution arising 
as an emergent phenomena. Furthermore, by allowing for arbitrary population 
sizes, FEArS can model population extinction. Arbitrary population sizes allows
us to simulate how a disease population within a patient may respond to 
therapy. In addition, FEArS models specific genotypes, allowing for more fine-
grained prediction of evolution.




# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

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

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References