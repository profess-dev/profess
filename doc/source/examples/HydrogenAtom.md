---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Hydrogen Atom

The hydrogen atom potential is $v(\mathbf{r})=-\tfrac{1}{r}$ in atomic units. Solving the associated Schrödinger equation yields the ground state energy $E_0=-\tfrac{1}{2}$.

For a single, isolated electron, we may infer its wave function from the electron density, $\psi(\mathbf{r}) = \sqrt{n(\mathbf{r})}$. Furthermore, its kinetic energy is given *exactly* by the Weizsäcker functional,

$$
T_W[n] = \int \mathrm{d}\mathbf{r} \, \sqrt{n(\mathbf{r})} \left(-\frac{1}{2}\nabla^2\right) \sqrt{n(\mathbf{r})}.
$$

For this reason, we can obtain $E_0$ as a simple example of orbital-free density functional theory,

$$
E_0 = \mathrm{min}_{n} \left[ T_{W}[n] + \int \mathrm{d}\mathbf{r} \, v(\mathbf{r}) n(\mathbf{r}) \right]
\quad\text{with}\quad
\int \mathrm{d}\mathbf{r} \, n(\mathbf{r}) = 1.
$$

The following code implements this approach, showing that the energy indeed approaches $-\frac{1}{2}$ when the simulation box is sufficiently large.

```{code-cell} ipython3
import profess
import matplotlib.pyplot as plt
import numpy as np

planewave_cutoff_energy = 50
box_lengths = np.linspace(2,20,10)

energies = []
for box_length in box_lengths:
    system = (
        profess.System.create_system(box_length*np.identity(3), planewave_cutoff_energy, ['b','h'])
        .add_coulomb_ions(1.0, [[0,0,0]], cutoff=0.5*box_length)
        .add_electrons()
        .add_weizsaecker_functional()
        .add_ion_electron_functional()
    )
    system.minimize_energy()
    energies.append(system.energy('h'))

plt.plot(box_lengths, energies)
plt.plot(box_lengths, -0.5*np.ones(box_lengths.size), 'k--')
plt.title('Hydrogen Atom Energy')
plt.xlabel('Box Length [Bohr]')
plt.ylabel('Energy [Hartree]')
plt.show()
```

One technical point: rather than the full Coulomb potential, the code uses a truncated version set to zero beyond a cutoff radius. This modification facilitates study of the isolated atom with periodic boundary conditions.
