---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  name: python3
---
# User Guide

```{toctree}
:maxdepth: 1
:hidden:
user-guide/build
user-guide/modify
user-guide/compute
```

Most users will only requre the ``profess.System`` class. The relevant commands fall into three categories:

* {doc}`building a new system <user-guide/build>`;
* {doc}`modifying an existing system <user-guide/modify>`;
* {doc}`computing or reporting system properties <user-guide/compute>`.


## Basic example

For a simple fcc aluminum crystal, the build stage might appear as follows.

```{code-cell} ipython3
:tags: [remove-input]
import matplotlib.pyplot as plt
import numpy as np

import pydeft as deft
import profess
```

```{code-cell} ipython3
:tags: [remove-output]

# choose box vectors and a planewave energy cutoff
box_vectors = 4.05 * np.identity(3)
energy_cutoff = 1200

# create the system
shape = profess.System.get_shape(box_vectors, energy_cutoff, ('a','ev'))
system = profess.System(shape)
system.set_box(box_vectors, 'a')

# add ions and electrons
system.add_ions(
    'examples/al.gga.recpot',
    box_vectors[0,0] * np.array([(0.0,0.0,0.0),
                                 (0.5,0.5,0.0),
                                 (0.5,0.0,0.5),
                                 (0.0,0.5,0.5)]),
    'a')
system.add_electrons(system.total_ion_charge())

# add the energy functional
(
system
    .add_wang_teter_functional()
    .add_hartree_functional()
    .add_perdew_burke_ernzerhof_functional()
    .add_ion_electron_functional()
)
```

The usual next step is to modify the system with {py:meth}`~profess.System.minimize_energy`,
finding the ground state electron density and energy.

```{code-cell} ipython3
:tags: [remove-output]
system.minimize_energy()
```

One can then print system properties, inspect the electron density, and more.

```{code-cell} ipython3
print('Volume (A3):        {:5.2f}' .format(system.volume('a3')))
print('Grid shape:         {:}'     .format(system.grid_shape))
print('Total energy (eV):  {:5.3f}' .format(system.energy('ev')))
```

```{code-cell}
fig, ax = plt.subplots(3, 1, constrained_layout=True, sharey=True)
fig.suptitle('Electron Density')
ax[0].plot([system.electron_density[i,0,0] for i in range(system.grid_shape[0])])
ax[1].plot([system.electron_density[i,i,0] for i in range(system.grid_shape[0])])
ax[2].plot([system.electron_density[i,i,i] for i in range(system.grid_shape[0])])
for i, label in enumerate(('[100]', '[110]', '[111]')):
    ax[i].set_xlabel(label + ' direction')
    ax[i].xaxis.set_ticklabels([])
    ax[i].set_ylabel('[a.u.]')
plt.show()
```

## Next steps

To see the full range of commands and capabilities, follow the links in the left sidebar.

Alternatively, browse the {doc}`examples`.
