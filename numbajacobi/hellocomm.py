import numba, numba_mpi, numpy

@numba.njit()
def hello():
  print(f"Rank: {numba_mpi.rank()}\nSize: {numba_mpi.size()}")

  src = numpy.array([1., 2., 3., 4., 5.])
  dst_tst = numpy.empty_like(src)

  if numba_mpi.rank() == 0:
    numba_mpi.send(src, dest=1, tag=11)
  elif numba_mpi.rank() == 1:
    numba_mpi.recv(dst_tst, source=0, tag=11)
  print(dst_tst)
hello()
