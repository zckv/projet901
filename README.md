# Lattice Bolltzam stencil with mpi4py and numba

## Command

mpiexec -n 2 poetry run python pylattice/rewrite.py pylattice/input_128x128.params files/obstacles_128x128.dat
python files/check/check.py --ref-av-vels-file=files/av_vels.dat --ref-final-state-file=files/final_state.dat --av-vels-file=av_vels.dat --final-state-file=final_state.dat

## LINKS

https://github.com/whdlgp/MPI_jacobi_iteration_example/blob/b80f4a66f5607df90f71aeeadef220f4769c4dd7/jacobi_mpi.py
https://github.com/mpi4py/mpi4py/blob/master/demo/mpi-ref-v1/ex-2.32.py
https://pypi.org/project/numba-mpi/
https://github.com/atmos-cloud-sim-uj/numba-mpi/tree/main/numba_mpi

