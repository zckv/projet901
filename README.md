# Lattice Bolltzam stencil with mpi4py and numba

## Installation

MPI est néccéssaire pour l'installation.

Avec poetry:
```
poerty install -n
```

Il est possible d'installer le projet avec les modules pip et venv.


## Lancement du programme

Pour lancer le programme, il y a deux possibilitées:

Avec poetry, en séquentiel:
```
poetry run python pylattice/d2q9-bgk.py pylattice/input_128x128.params pylattice/obstacles_128x128.dat
```


Avec mpiexec:
```
mpiexec -n X poetry run python pylattice/rewrite.py pylattice/input_128x128.params pylattice/obstacles_128x128.dat
```
X le nombres de processus à lancer

Il est possible de vérifier les résultats avec le script fourni:
```
poetry run python files/check/check.py --ref-av-vels-file=files/av_vels.dat --ref-final-state-file=files/final_state.dat --av-vels-file=av_vels.dat --final-state-file=final_state.dat
```

