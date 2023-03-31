from numba import njit, stencil, prange, float64, int32, boolean
from argparse import ArgumentParser
from math import sqrt, ceil
from sys import argv

import coloredlogs
import numba_mpi as nm
import logging
import numpy as np
import time
import os
import sys

NSPEEDS = 9
FINALSTATEFILE = "final_state.dat"
INITIALSTATEFILE = "initial_state.dat"
AVVELSFILE = "av_vels.dat"

llevel = logging.INFO

logging.basicConfig(encoding='utf-8', level=llevel)
logger = logging.getLogger(os.path.basename(sys.argv[0]))
logger.setLevel(llevel)
coloredlogs.DEFAULT_FIELD_STYLES["levelname"]["color"] = "cyan"
coloredlogs.install(logger=logger, level=llevel)


def read_input_file(file_path):
    d = {}
    with open(file_path, 'r') as f:
        d["nx"] = int(f.readline())
        d["ny"] = int(f.readline())
        d["maxIters"] = int(f.readline())
        d["reynolds_dim"] = int(f.readline())
        d["density"] = float(f.readline())
        d["accel"] = float(f.readline())
        d["omega"] = float(f.readline())
    return d


def read_obstacle_file(file_path):
    block = set()
    with open(file_path, 'r') as f:
        for line in f:
            x, y, b = line.split()
            block.add((int(x), int(y)))
    return block


def main():
    # TODO initalization
    param = read_input_file(argv[1])
    block = read_obstacle_file(argv[2])

    tot_tic = time.time()
    init_tic = time.time()

    w0 = param["density"] * 4 / 9
    w1 = param["density"] / 9
    w2 = param["density"] / 36

    nx = param["nx"]
    ny = param["ny"]
    density = param["density"]
    accel = param["accel"]
    omega = param["omega"]
    reynolds_dim = param["reynolds_dim"]

    # main_grid
    cells = np.array([[[w0, w1, w1, w1, w1, w2, w2, w2, w2] for _ in range(nx)] for _ in range(ny)], dtype=np.float64)
    obstacles = np.zeros((ny, nx), dtype=np.bool_)
    for x, y in block:
            obstacles[x,y] = 1
    logger.debug(param)

    init_toc = time.time()
    comp_tic = time.time()
    # logger.info(f"Rank: {nm.rank()}; {start_blk} {end_blk}")

    av_vels = compute(cells, obstacles, param["maxIters"], nx, ny, density, accel, omega)

    comp_toc = time.time()
    col_tic = time.time()

    collat_cells(cells, nx)

    col_toc = time.time()
    tot_toc = time.time()

    if nm.rank() == 0:
        print(f"Reynolds number:\t\t{calc_reynolds(cells, obstacles, nx, ny, omega, reynolds_dim)}")
        print(f"Elapsed Init time:\t\t{init_toc - init_tic} seconds")
        print(f"Elapsed Compute time:\t\t{comp_toc - comp_tic} seconds")
        print(f"Elapsed Collate time:\t\t{col_toc - col_tic} seconds")
        print(f"Elapsed Total time:\t\t{tot_toc - tot_tic} seconds")

        write_values(cells, obstacles, av_vels, nx, ny, density, param["maxIters"])



@njit(float64(float64), parallel=False)
def exchange_av_vel(av_vel_t):
    av_vels = np.zeros((nm.size()), dtype=np.float64)
    av_vels[nm.rank()] = av_vel_t
    if nm.size() > 1:
        if nm.rank()==0:
            buff = np.zeros((nm.size()), dtype=np.float64)
            for j in range(1, nm.size()):
                nm.recv(buff, source=j, tag=40)
                av_vels[j] = buff[j]
        else:
            nm.send(av_vels, dest=0, tag=40)
    return av_vels.mean()


@njit((float64[:, :, :], int32), parallel=False)
def exchange_halos(cells, nx):
    blk_sz = nx//nm.size()
    if nm.size() > 1:
        for i in range(nm.size()):
            start_blk = i * blk_sz
            end_blk = i * blk_sz + blk_sz - 1
            if i == nm.rank():
                for j in range(nm.size()):
                    if j != i:
                        nm.send(cells[start_blk], dest=j, tag=10)
                        nm.send(cells[end_blk], dest=j, tag=20)
            else:
                nm.recv(cells[start_blk], source=i, tag=10)
                nm.recv(cells[end_blk], source=i, tag=20)

@njit((float64[:, :, :], int32), parallel=False)
def collat_cells(cells, nx):
    blk_sz = nx//nm.size()
    if nm.size() > 1:
        start_blk = nm.rank() * blk_sz
        end_blk = nm.rank() * blk_sz + blk_sz -1
        if nm.rank()==0:
            for j in range(1, nm.size()):
                j_st = j*blk_sz
                j_end = j*blk_sz + blk_sz - 1
                nm.recv(cells[j_st:j_end+1], source=j, tag=30)
        else:
            nm.send(cells[start_blk:end_blk + 1], dest=0, tag=30)


@njit(int32(float64[:, :, :], int32, int32), parallel=True)
def total_density(cells, nx, ny):
    return cells.sum()


@njit(float64(float64[:,:,:], boolean[:,:], int32, int32, int32, int32), parallel=True)
def av_velocity(cells, obstacles, nx, ny, start_blk, end_blk,):
    """TODO"""
    tot_cells = 0
    tot_u = 0.

    for jj in range(ny):
        for ii in range(start_blk, end_blk+1):
            if not obstacles[ii][jj]:
                c = cells[ii][jj]
                ld =  c.sum()
                u_x = (c[1] + c[5] + c[8] - (c[3] + c[6] + c[7])) / ld
                u_y = (c[2] + c[5] + c[6] - (c[4] + c[7] + c[8])) / ld
                tot_u += np.sqrt((u_x * u_x) + (u_y * u_y))
                tot_cells += 1
    return tot_u / tot_cells


@njit(float64(float64[:, :, :], boolean[:, :], int32, int32, float64, int32), parallel=False)
def calc_reynolds(cells, obstacles, nx, ny, omega, reynolds_dim):
    viscosity = 1. / 6. * (2. / omega - 1.)
    return av_velocity(cells, obstacles, nx, ny, 0, nx-1) * reynolds_dim / viscosity


@njit((float64[:,:,:], boolean[:,:], int32, int32, int32, int32, float64, float64, float64), parallel=True)
def timestep(
    cells,
    obstacles,
    start_blk,
    end_blk,
    nx,
    ny,
    density,
    accel,
    omega
):
    """One step elapse"""
    c_sq = 1. / 3.  # square of speed of sound
    w0 = 4. / 9.  # weighting factor
    w1 = 1. / 9.  # weighting factor
    w2 = 1. / 36.  # weighting factor
    w3 = density * accel / 9.
    w4 = density * accel / 36.

    tmp = cells.copy()

    for jj in prange(ny):
        for ii in prange(start_blk, end_blk + 1):
            # init in parallel sec
            c = np.zeros((NSPEEDS), dtype=np.float64)
            u = np.zeros((NSPEEDS), dtype=np.float64)

            # accelerate and propagate
            y_n = (jj + 1) % ny
            x_e = (ii + 1) % nx
            y_s = (ny - 1) if jj == 0 else jj - 1 #  * ceil((jj)/(jj+1)) +  (jj - 1) * ceil((jj)/(jj+1) + 1 % 2)
            x_w = (nx - 1) if ii == 0 else ii - 1

            w3c = (w3 if jj == ny -2 else 0)
            w4cn = (w4 if y_n == ny - 2 else 0)
            w4cs = (w4 if y_s == ny - 2 else 0)

            c[0] = tmp[ii][jj][0]
            c[1] = tmp[x_e][jj][1] + w3c
            c[2] = tmp[ii][y_n][2]
            c[3] = tmp[x_w][jj][3] - w3c
            c[4] = tmp[ii][y_s][4]
            c[5] = tmp[x_e][y_n][5] + w4cn
            c[6] = tmp[x_w][y_n][6] - w4cn
            c[7] = tmp[x_w][y_s][7] - w4cs
            c[8] = tmp[x_e][y_s][8] + w4cs

            # rebound and collision
            if obstacles[ii, jj]:
                cells[ii][jj][1] = c[3]
                cells[ii][jj][2] = c[4]
                cells[ii][jj][3] = c[1]
                cells[ii][jj][4] = c[2]
                cells[ii][jj][5] = c[7]
                cells[ii][jj][6] = c[8]
                cells[ii][jj][7] = c[5]
                cells[ii][jj][8] = c[6]
            else:
                ld =  c.sum()
                u_x = (c[1] + c[5] + c[8] - (c[3] + c[6] + c[7])) / ld
                u_y = (c[2] + c[5] + c[6] - (c[4] + c[7] + c[8])) / ld
                u_sq = u_x * u_x + u_y * u_y
                u_sq22 = 2 * c_sq * c_sq
                uc_sq = u_sq / (2*c_sq)

                u[1] = u_x / c_sq + (u_x * u_x) / u_sq22
                u[2] = u_y / c_sq + (u_y * u_y) / u_sq22
                u[3] = -u_x / c_sq + (-u_x * -u_x) / u_sq22
                u[4] = -u_y / c_sq + (-u_y * -u_y) / u_sq22
                u[5] = (u_x + u_y) / c_sq + ((u_x + u_y) * (u_x + u_y)) / u_sq22
                u[6] = (-u_x + u_y) / c_sq + ((-u_x + u_y) * (-u_x + u_y)) / u_sq22
                u[7] = (-u_x - u_y) / c_sq + ((-u_x - u_y) * (-u_x - u_y)) / u_sq22
                u[8] = (u_x - u_y) / c_sq + ((u_x - u_y) * (u_x - u_y)) / u_sq22

                # equilibrium densities
                # zero velocity density: weight w0
                cells[ii][jj][0] = c[0] + omega * (w0 * ld * (1. - uc_sq) - c[0])

                # axis speeds: weight w1 */
                cells[ii][jj][1] = c[1] + omega * (w1 * ld * (1 + u[1] - uc_sq) - c[1])
                cells[ii][jj][2] = c[2] + omega * (w1 * ld * (1 + u[2] - uc_sq) - c[2])
                cells[ii][jj][3] = c[3] + omega * (w1 * ld * (1 + u[3] - uc_sq) - c[3])
                cells[ii][jj][4] = c[4] + omega * (w1 * ld * (1 + u[4] - uc_sq) - c[4])

                # diagonal speeds: weight w2 */
                cells[ii][jj][5] = c[5] + omega * (w2 * ld * (1 + u[5] - uc_sq) - c[5])
                cells[ii][jj][6] = c[6] + omega * (w2 * ld * (1 + u[6] - uc_sq) - c[6])
                cells[ii][jj][7] = c[7] + omega * (w2 * ld * (1 + u[7] - uc_sq) - c[7])
                cells[ii][jj][8] = c[8] + omega * (w2 * ld * (1 + u[8] - uc_sq) - c[8])


@njit(float64[:](float64[:,:,:], boolean[:,:], int32, int32, int32, float64, float64, float64), parallel=False)
def compute(cells, obstacles, n, nx, ny, density, accel, omega):
    av_vels = np.zeros(n, dtype=np.float64)
    blk_sz = nx//nm.size()
    start_blk = np.int32(nm.rank()*blk_sz)
    end_blk = np.int32(nm.rank()*blk_sz + blk_sz - 1)

    for tt in range(n):
        timestep(
            cells,
            obstacles,
            start_blk,
            end_blk,
            nx,
            ny,
            density,
            accel,
            omega
        )
        exchange_halos(cells, nx)
        av_vels[tt] = av_velocity(cells, obstacles, nx, ny, start_blk, end_blk)
        av_vels[tt] = exchange_av_vel(av_vels[tt])
        # if nm.rank() == 0:
        #    logger.info(f"==timestep: {tt}==")
        #    logger.info(f"ev velocity: {av_vels[tt]}")
        #    logger.info(f"tot density: {total_density(cells, nx, ny)}")
    return av_vels



def write_values(cells, obstacles, av_vels, nx, ny, density, maxIters):
    c_sq = 1. / 3.
    with open(FINALSTATEFILE, "w") as f:
        for jj in range(ny):
            for ii in range(nx):
                if obstacles[ii, jj]:
                    u_x = 0
                    u_y = 0
                    u = 0

                    pressure = density * c_sq
                else:
                    local_density = cells[ii][jj].sum()
                    # compute x velocity component
                    u_x = (
                                 cells[ii][jj][1] + cells[ii][jj][5] +
                                 cells[ii][jj][8] - (
                                         cells[ii][jj][3] +
                                         cells[ii][jj][6] +
                                         cells[ii][jj][7]
                                 )
                         ) / local_density
                    # compute y velocity component
                    u_y = (
                                  cells[ii][jj][2] +
                                  cells[ii][jj][5] +
                                  cells[ii][jj][6] - (
                                          cells[ii][jj][4] +
                                          cells[ii][jj][7] +
                                          cells[ii][jj][8]
                                  )
                          ) / local_density

                    # compute norm of velocity
                    u = sqrt((u_x * u_x) + (u_y * u_y))
                    # compute pressure
                    pressure = local_density * c_sq
                f.write(f"{ii} {jj} {u_x} {u_y} {u} {pressure} {int(obstacles[ii, jj])}\n")

    with open(AVVELSFILE, "w") as f:
        for ii in range(maxIters):
            f.write(f"{ii} {av_vels[ii]}\n")



if __name__ == "__main__":
    main()
