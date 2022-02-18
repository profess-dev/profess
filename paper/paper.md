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
  - name: Department of Mechanical and Aerospace Engineering, Princeton University
    index: 4
  - name: Department of Chemical and Biomolecular Engineering, University of California, Los Angeles
    index: 5
  - name: Office of the Chancellor, University of California, Los Angeles
    index: 6
date: 22 August 2021
bibliography: paper.bib
---

# Summary

A core endeavor of computational physics is to determine the properties of materials from quantum-mechanical first principles. Density functional theory (DFT) has achieved remarkable success in this effort, providing both a rigorous framework and pragmatic approximations applicable across material classes. The `PROFESS` code (“PRinceton Orbital-Free Electronic Structure Software”) implements orbital-free DFT, which uses the electron density alone—bypassing wave functions and Schrödinger's equation entirely—to characterize a many-electron system. This formulation, while imposing a few extra challenges, remains exact in principle and offers an especially elegant and computationally efficient solution to the electronic structure problem.

The new version of `PROFESS` described here, completely rewritten, is a high-performance engine for orbital-free DFT that will promote its rapid development and adoption.

# Statement of need

The materials modelling software ecosystem confronts several disparate aims. Its simulations frequently operate at computing frontiers, demanding optimized methods and code. At the same time, it must incorporate theoretical innovations emerging from many subfields at once, creating strong preferences for modularity and interoperability, as well as a need for efficient prototyping workflows. While `PROFESS` [@ho_introducing_2008; @hung_introducing_2010;  @chen_introducing_2015] has achieved demonstrable success for research applications in materials science [@witt_orbital-free_2018], it was conceived as a self-contained entity. Prior versions lacked agility, constrained by an older paradigm.

This fourth version of `PROFESS` is completely new, although it draws on the innovations of earlier releases. Its main features are implemented in modern C++, compiled as a library. Deep access from both Python and Julia is enabled by `pybind11` and `Cxxwrap.jl`, respectively. This structure, without sacrificing performance, facilitates rapid scripting of common tasks, as well as frictionless partnerships with complementary tools in the materials modelling software ecosystem.

* _Interoperability_. Over the past decade, scripting languages like Python have emerged as essential "glue languages," particularly for materials modelling. For example, the Atomic Simulation Environment and Phonopy are two Python-based tools that, when used in conjunction with a DFT code, provide robust geometry optimizations and many other features that are common requirements for DFT codes. Historically, `PROFESS` has been used. A chief example is the `PROFESS@QE` interface, which enables. At the same time, `PROFESS` has been used in conjunction with Phonopy and ATAT. However, older versions of `PROFESS` were predominantly Fortran code, compiled into a standalone executable, these interfaces were frequently brittle and carried an unwelcome maintence burden. We expect that the new design, particularly the clean API and scripting mode, will promote new partnerships of this kind, and will lower the maintence burden.

* _Unit testing_. A second clear benefit of the new structure is ease of unit testing. Because the C++ classes are individually wrapped and accessible from Python, it is simple to create granular tests of individual components. Moreover, such tests can form an overall framework that facilitates continuous integration.

* _Rapid prototyping_. Frequently, we have found it easier to develop and test new features in Python, before ultimately translating them into compiled code. Such hybrid workflows were generally impossible with older versions of `PROFESS`.


> But, its full potential is constrained by its design;
conceived as a self-contained entity, it is not readily interoperable with other CMS tools. An example illustrates the point: several years ago, another research group (not close collaborators) produced a PROFESS@QuantumEspresso hybrid59 to study the hard-to-access warm dense matter regime with OFDFT functionals—embodying precisely the creative fusion advocated for here. However, the task required an unwelcome amount of hacking, and changes to the progenitor codes have not necessarily propagated gracefully to the hybrid. A central part of the MolSSI mission is to remove impediments of this
kind, which can “[delay] by years the practical realization of theoretical innovations.”60

> A key task, then—the core of the proposal—is to eliminate any barriers that would impede this researcher’s access to OFDFT energies, forces, and other properties. PROFESS computes first-principles energies (as well as forces and stresses), conducts structural relaxations, and performs MD simulations. 

> completely new
main codebase is modern C++
deep access to the C++ classes is provided with aid from the pybind11 wrapping tool

>This strategy confers a few benefits:




# Background

As authorized by the foundational theorems of density functional theory, orbital-free DFT finds the ground state energy, $E_0$, by solving a minimization problem,
\begin{equation} \label{eq:E0}
E_0 = \min_n E[n].
\end{equation}
The function $n(\mathbf{r})$ represents the electron density and $E[n]$ is a functional that returns an energy value all candidate densities. The minimization in \autoref{eq:E0}, which yields the ground state energy, $E_0$, is performed over all admissible density functions, which must have $n(\mathbf{r}) \ge 0$ everywhere and which must integrate to $N$, the total number of electrons.
One simple method of enforcing such constraints are to write the density as $n(\mathbf{r}) = (N/\bar{N}) \chi^2(\mathbf{r})$, where $\chi$ is an unconstrained auxiliary function and $\bar{N}=\int \chi^2(\mathbf{r}') \, d\mathbf{r}'$. The functional derivative of $E[n]$ with respect to $\chi$, is then
\begin{equation} \label{eq:dE_dchi}
\frac{\delta E}{\delta \chi(\mathbf{r})} = \frac{N}{\bar{N}} \left[v(\mathbf{r}) - \frac{1}{N} \int v(\mathbf{r}) n(\mathbf{r}) \right] ,
\end{equation}
where $v(\mathbf{r})$ is the potential associated with the energy functional: $v(\mathbf{r}) = \delta E / \delta n(\mathbf{r})$. \autoref{eq:dE_dchi} facilitates minimization of the energy with respect to $\chi(\mathbf{r}$.

# Example application

mention random structure searching too


# detritus

For an ever-increasing set of materials, problems involving length and time scales inaccessible with competing methods can become routine. 

* Completely rewritten
* Primarily C++ with Python wrappers for everything
* Prioritizing modularity
    *  DEFT for math operations
    *  Separate things (like Ewald sum) teased out
    *  Functionals easily loaded into Python
    * libKEDF
* Prioritizing interoperability
    *  Many examples
* GPU


# Acknowledgements

Wojciech Jankowski, Chuin Wei Tan, Thomas Ginnis, Pascal Salzbrenner
Johannes Dieterich, Florian Libisch
MolSSI

# References
