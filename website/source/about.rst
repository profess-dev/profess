*****
About
*****

``profess`` began as PRinceton Orbital-Free Electronic Structure Software in
the research group of `Prof. Emily Carter
<https://research.seas.ucla.edu/carter/>`_. Over time, the code and
contributors grew beyond Princeton, and ``profess`` is now an open-source
project independent of any single institution.

The current version of ``profess`` includes a C++ library of core features, the
``profess`` executable for launching simple calculations, and Python bindings
that enable more complex workflows. This architecture ensures high performance,
while enhancing interoperability with other tools. Older versions of
``profess`` (including ``profess`` 3, available `here
<https://github.com/EACcodes/PROFESS>`_) were mostly written in Fortran.

History
-------

.. sidebar:: Linear-scaling parallel algorithms for the first principles treatment of metals

    S\. C. Watson and E. A. Carter, Computer Physics Communications **128**, 67 (2000)

    DOI: `10.1016/S0010-4655(00)00064-3 <https://doi.org/10.1016/S0010-4655(00)00064-3>`_

The Carter Group's work on orbital-free density functional theory (OFDFT) began
ca. 1996, coinciding with a sabbatical visit to Paul Madden's group at Oxford.
Afterwards, Stuart Watson, having completed his PhD in the Madden Group, joined
the Carter Group in UCLA as a postdoc. This collaboration produced a paper
summarizing the key algorithms for a planewave OFDFT code.

The first official ``profess`` release was published in 2008, after the Carter
Group's move to Princeton, with principal authors Greg Ho and Vincent Lignères.

.. topic:: Introducing PROFESS: A new program for orbital-free density functional theory calculations

    G\. S. Ho, V. L. Lignères, and E. A. Carter, Computer Physics Communications **179**, 839 (2008)

    DOI: `10.1016/j.cpc.2008.07.002 <https://doi.org/10.1016/j.cpc.2008.07.002>`_

.. sidebar:: Accurate simulations of metals at the mesoscale: Explicit treatment of 1 million atoms with quantum mechanics

    L\. Hung and E. A. Carter, Chemical Physics Letters 475, 163 (2009)

    DOI: `10.1016/j.cplett.2009.04.059 <https://doi.org/10.1016/j.cplett.2009.04.059>`_

Subsequent work by Linda Hung enabled the major milestone of million-atom
calculations. It also led to a second ``profess`` release, which added Chen
Huang and Issac Shin as new contributors.

.. topic:: Introducing PROFESS 2.0: A parallelized, fully linear scaling program for orbital-free density functional theory calculations

    L\. Hung, C. Huang, I. Shin, G. S. Ho, V. L. Lignères, and E. A. Carter, Computer Physics Communications **181**, 2208 (2010)

    DOI: `10.1016/j.cpc.2010.09.001 <https://doi.org/10.1016/j.cpc.2010.09.001>`_

Further advances by Mohan Chen, Junchao Xia, and Johannes Dieterich, together
with others already mentioned, gave rise to ``profess`` 3, featuring molecular
dynamics with OFDFT.

.. topic:: Introducing PROFESS 3.0: An advanced program for orbital-free density functional theory molecular dynamics simulations

    M\. Chen, J. Xia, C. Huang, J. M. Dieterich, L. Hung, I. Shin, and E. A. Carter, Computer Physics Communications **190**, 228 (2015)

    DOI: `10.1016/j.cpc.2014.12.021 <https://doi.org/10.1016/j.cpc.2014.12.021>`_

.. sidebar:: libKEDF: An accelerated library of kinetic energy density functionals

    J\. M. Dieterich, W. C. Witt, and E. A. Carter, Journal of Computational Chemistry **38**, 1552 (2017)

    DOI: `10.1002/jcc.24806 <https://doi.org/10.1002/jcc.24806>`_

The transition to C++ began with work by Johannes Dieterich and Chuck Witt, who
produced ``libKEDF``, a library of nonlocal kinetic energy functionals. The
success of that experiment was one motivation for the complete redesign of the
code base. The current release, ``profess 4``, prioritizes interoperability
with other materials modelling tools.
