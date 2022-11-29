---
title: '`PROFESS` 4: New Horizons for Orbital-free Density Functional Theory'
tags:
  - computational physics
  - density functional theory
  - orbital-free density functional theory
authors:
  - name: William C. Witt^[corresponding author]
    orcid: 0000-0002-1578-1888
    affiliation: "1, 2"
  - name: Chris J. Pickard
    orcid: 0000-0002-9684-5432
    affiliation: "1, 3"
  - name: Emily A. Carter
    orcid: 0000-0001-7330-7554
    affiliation: "4, 5, 6"
affiliations:
  - name: Department of Materials Science and Metallurgy, University of Cambridge
    index: 1
  - name: Christ's College, University of Cambridge
    index: 2
  - name: Advanced Institute for Materials Research, Tohoku University
    index: 3
  - name: Department of Mechanical and Aerospace Engineering, Andlinger Center for Energy and the Environment, and Department of Applied and Computational Mathematics, Princeton University
    index: 4
  - name: Princeton Plasma Physics Laboratory
    index: 5
date: 22 August 2021
bibliography: paper.bib
---

# Summary

A core endeavor of computational physics is to determine the properties of materials from quantum-mechanical first principles. Density functional theory (DFT) has achieved remarkable success in this effort, providing both a rigorous framework and pragmatic approximations applicable across material classes. The `PROFESS` code (“PRinceton Orbital-Free Electronic Structure Software”) implements orbital-free DFT, which uses the electron density alone—bypassing Schrödinger's equation and its associated wave functions entirely—to characterize a many-electron system. This formulation, while imposing a few extra challenges, remains exact in principle and offers an especially elegant and computationally efficient solution to the electronic structure problem.

The new version of `PROFESS`, completely rewritten, is a high-performance engine for orbital-free DFT that will facilitate its rapid development and adoption.

# Background

Using the foundational theorems of density functional theory (DFT) [@hohenberg_inhomogeneous_1964; @kohn_self-consistent_1965], `PROFESS` finds the ground state energy, $E_0$, of a many-electron system by solving a minimization problem,
\begin{equation} \label{eq:E0}
E_0 = \min_n E[n].
\end{equation}
The function $n(\mathbf{r})$ represents the electron density and $E[n]$ is a functional that returns an energy value for all candidate densities, obtaining the total energy as a sum of kinetic, electrostatic, and exchange-correlation contributions.

The minimization in \autoref{eq:E0} is performed over all admissible density functions, which must have $n(\mathbf{r}) \ge 0$ everywhere and which must integrate to $N$, the total number of electrons. One method of enforcing these constraints (used by `PROFESS` 4) is to write the density as $n(\mathbf{r}) = (N/\bar{N}) \chi^2(\mathbf{r})$, where $\chi$ is an unconstrained auxiliary function and $\bar{N}=\int \chi^2(\mathbf{r}') \, d\mathbf{r}'$. The functional derivative of $E[n[\chi]]$ with respect to $\chi$, is then
\begin{equation} \label{eq:dE_dchi}
\frac{\delta E}{\delta \chi(\mathbf{r})} = \frac{N}{\bar{N}} \left[\frac{\delta E}{\delta n(\mathbf{r})} - \underbrace{\frac{1}{N} \int v(\mathbf{r}) n(\mathbf{r})}_{\mu} \right] .
\end{equation}
Using \autoref{eq:dE_dchi}, it is straightforward to minimize the energy with, for example, a quasi-Newton method, after which one may identify the rightmost term in the brackets as the Fermi energy, $\mu$.

# Statement of need

Materials modelling software must achieve several disparate aims. Its simulations frequently operate at computing frontiers, demanding optimized methods and code. At the same time, it must incorporate theoretical innovations emerging from many subfields at once, creating strong preferences for modularity and interoperability, as well as a need for efficient prototyping workflows. While `PROFESS` [@ho_introducing_2008; @hung_introducing_2010;  @chen_introducing_2015] has achieved demonstrable success for research applications in materials science [@witt_orbital-free_2018], it was conceived as a self-contained entity. Prior versions lacked agility, constrained by an older paradigm.

This fourth version of `PROFESS` is completely new, although it draws on the innovations of earlier releases. Its main features are implemented in modern C++ code, compiled into a library. Deep access from both Python and Julia is enabled by `pybind11` and `CxxWrap.jl`, respectively. This structure, without sacrificing performance, facilitates rapid scripting of common tasks, as well as frictionless partnerships with complementary tools. Benefits include:

* _Interoperability_. Over the past decade, Python and (more recently) Julia have emerged as essential "glue languages." For example, the Atomic Simulation Environment [@larsen_atomic_2017] and Phonopy [@togo_first_2015] are two widely used Python-based tools that, when coupled to a DFT code, provide robust geometry optimizations and many other features. As older versions of `PROFESS` were predominantly Fortran code compiled into a standalone executable, interfaces to such tools were generally brittle and carried an unwelcome maintenance burden. The new design, particularly the library interface and scripting modes, eliminates these barriers.

* _Unit testing_. A second clear benefit of the new structure is ease of unit testing. Because the C++ classes are individually wrapped and accessible from Python, it is simple to create granular tests of individual components. Moreover, such tests can form an overall framework that facilitates continuous integration.

* _Rapid prototyping_. It is often simpler (and less error-prone) to develop and test extensions in a scripting language, before ultimately translating them into compiled code. Such hybrid workflows were generally impossible with older versions of `PROFESS`, but are now routine.

All together, the new `PROFESS` structure improves usability and facilitates rapid development, while retaining the advantages of compiled code. For example, recent work using orbital-free DFT for crystal structure prediction was enabled by the robust geometry optimizations offered by the Atomic Simulation Environment and would not have been possible with earlier `PROFESS` releases.

# Acknowledgements

The authors gratefully acknowledge Wojciech Jankowski, Chuin Wei Tan, Thomas Ginnis, Pascal Salzbrenner, Johannes Dieterich, and Florian Libisch.

# References
