DFT output data used for producing the figures.

All structures and energies are provided in the ```.traj```-format of the atomic simulation environment (ASE).

Note that all energies correspond to constant potential energies resulting from constant potential calculations within the [Solvated Jellium method](https://wiki.fysik.dtu.dk/gpaw/documentation/sjm/sjm.html)

Both the number of electrons and potential are included in ```results``` dictionary within the trajectory files.

WARNING: Sadly, in order to parse the trajectory files properly a slightly changed ASE version from trunk is needed, as vanilla ASE does not allow readinig and writing of the number of electrons and potentials in the trajectory format. The necessary changes to ASE are:

In ase/calculators/singlepoint.py change line 25 from
```
if property in ['energy', 'magmom', 'free_energy']:
```
to
```
if property in ['energy', 'magmom', 'free_energy','ne','electrode_potential']:
```


and in ase/calculators/calculator line 98 exchange:
```
all_properties = ['energy', 'forces', 'stress', 'stresses', 'dipole',
                  'charges', 'magmom', 'magmoms', 'free_energy']
```
with
```
all_properties = ['energy', 'forces', 'stress', 'stresses', 'dipole',
                  'charges', 'magmom', 'magmoms', 'free_energy','ne','electrode_potential']
```

We also note that GPAW version 21.1.1b1 was used for the calculations and more recent version changed the names of `ne` and `electrode_potential`. However, ASE can still read the results if above's changes are made.

