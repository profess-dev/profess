import pydeft
import profess
import numpy as np

energy_cutoff = 1200
box_vectors = 4.050 * np.array([[0.5, 0.5, 0.0],
                                [0.0, 0.5, 0.5],
                                [0.5, 0.0, 0.5]])
system = profess.System.create(box_vectors, energy_cutoff, ['a','ev'])
system.add_ions('potentials/al.gga.recpot', np.zeros([1,3]), 'a')
system.add_ion_electron_functional()
print(system.external_potential())
