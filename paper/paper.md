---
title: 'PROFESS 4'
tags:
  - computational materials science
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
  - name: Department of Materials Science and Metallurgy, University of Cambridge, Cambridge, United Kingdom
    index: 1
  - name: Christ's College, University of Cambridge, Cambridge, United Kingdom
    index: 2
  - name: Advanced Institute for Materials Research, Tohoku University, Sendai, Japan
    index: 3
  - name: Department of Mechanical and Aerospace Engineering, Princeton University, Princeton, New Jersey, United States
    index: 4
  - name: Department of Chemical and Biomolecular Engineering, University of California, Los Angeles, Los Angeles, California, United States
    index: 5
  - name: Office of the Chancellor, University of California, Los Angeles, Los Angeles, California, United States
    index: 6
date: 22 August 2021
bibliography: paper.bib
---

# Summary

A key endeavor of computational physics is to predict and interpret the properties of materials from quantum-mechanical first principles. Density functional theory (DFT) has achieved remarkable success at this task, providing both a rigorous framework and pragmatic approximate implementations. The PROFESS code (“PRinceton Orbital-Free Electronic Structure Software”) implements orbital-free DFT, which, while simultaneously older and less well-developed than conventional DFT, offers an especially elegant and computationally efficient approach to the electronic structure problem. Something about the type of problems OF-DFT can solve.

The radical difference between orbital-free DFT and conventional DFT is that the orbital-free variant uses the electron density alone—no wave functions—to determine the energy of a many-electron system. For further distinguishing features of orbital-free DFT, see some review papers.

PROFESS  is a high-performance engine for orbital-free DFT.

### The core task in orbital-free DFT

As authorized by the foundational theorems of density functional theory, orbital-free DFT finds the ground state energy, $E_0$, by solving a minimization problem,
$$
E_0 = \min_n E[n].
$$
The function $n(\mathbf{r})$ represents the electron density and $E[n]$ is a functional that returns an energy value all candidate densities. The minimization, which yields the ground state energy, $E_0$, is performed over all admissible density functions, which must have $n(\mathbf{r}) \ge 0$ everywhere and which must integrate to $N$, the total number of electrons.
One simple method of enforcing such constraints are to write the density as $n(\mathbf{r}) = (N/\bar{N}) \chi^2(\mathbf{r})$, where $\chi$ is an unconstrained auxiliary function and $\bar{N}=\int \chi^2(\mathbf{r}') \, d\mathbf{r}'$. The functional derivative of $E[n]$ with respect to $\chi$, is then
$$
\frac{\delta E}{\delta \chi(\mathbf{r})} = \frac{N}{\bar{N}} \left[v(\mathbf{r}) - \frac{1}{N} \int v(\mathbf{r}) n(\mathbf{r}) \right] ,
$$
where $v(\mathbf{r})$ is the potential associated with the energy functional: $v(\mathbf{r}) = \delta E / \delta n(\mathbf{r})$.

### PROFESS 4

completely new
main codebase is modern C++
deep access to the C++ classes is provided with aid from the pybind11 wrapping tool

This strategy confers a few benefits:
* interoperability with other tools
    * "glue languages"
    * PROFESS@QE
    * phonopy and ATAT.
* unit testing
* rapid development


The PROFESS 4 code is completely new, although it draws on the innovations of earlier PROFESS releases. Its main features are implemented in modern C++; however, deep access from Python is provided with the pybind11 wrapping tool. This structure, without sacrificing performance, facilitates rapid scripting of common tasks, as well as frictionless partnerships with complementary tools in the materials modelling software ecosystem. The same philosophy is adopted, for example, by the Psi4 code.

* _Interoperability_. Over the past decade, scripting languages like Python and Julia have emerged as essential "glue languages," particularly for materials modelling. For example, the Atomic Simulation Environment and Phonopy are two Python-based tools that, when used in conjunction with a DFT code, provide robust geometry optimizations and many other features that are common requirements for DFT codes. Historically, PROFESS has been used. A chief example is the PROFESS@QE interface, which enables. At the same time, PROFESS has been used in conjunction with Phonopy and ATAT. However, older versions of PROFESS were predominantly Fortran code, compiled into a standalone executable, these interfaces were frequently brittle and carried an unwelcome maintence burden. We expect that the new design, particularly the clean API and scripting mode, will promote new partnerships of this kind, and will lower the maintence burden.

* _Unit testing_. A second clear benefit of the new structure is ease of unit testing. Because the C++ classes are individually wrapped and accessible from Python, it is simple to create granular tests of individual components. Moreover, such tests can form an overall framework that facilitates continuous integration.

* _Rapid prototyping_. Frequently, we have found it easier to develop and test new features in Python, before ultimately translating them into compiled code. Such hybrid workflows were generally impossible with older versions of PROFESS.

### Example application

mention random structure searching too


### notes

Completely rewritten
Primarily C++ with Python wrappers for everything
Prioritizing modularity
  DEFT for math operations
  Separate things (like Ewald sum) teased out
  Functionals easily loaded into Python
  libKEDF
Prioritizing interoperability
  Many examples
GPU


### Acknowledgements

Wojciech Jankowski, Chuin Wei Tan, Thomas Ginnis, Pascal Salzbrenner
Johannes Dieterich, Florian Libisch
MolSSI
