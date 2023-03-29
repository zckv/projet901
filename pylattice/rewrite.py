from numba import njit, stencil, prange
from argparse import ArgumentParser
from math import sqrt
from sys import argv

import coloredlogs
import logging
import numpy as np
import time
import os
import sys

NSPEEDS = 9
FINALSTATEFILE = "final_state.dat"
INITIALSTATEFILE = "initial_state.dat"
AVVELSFILE = "av_vels.dat"

logging.basicConfig(encoding='utf-8', level=logging.DEBUG)
logger = logging.getLogger(os.path.basename(sys.argv[0]))
logger.setLevel(logging.INFO)
coloredlogs.DEFAULT_FIELD_STYLES["levelname"]["color"] = "cyan"
coloredlogs.install(logger=logger, level=logging.INFO)


class ProjectParser(ArgumentParser):
    """Parse args from command line"""
    # TODO


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

    logger.debug("==start program==")

    param = read_input_file(argv[1])
    block = read_obstacle_file(argv[2])

    tot_tic = time.time_ns()
    init_tic = time.time_ns()

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
    cells = np.array([[[w0, w1, w1, w1, w1, w2, w2, w2, w2] for _ in range(nx)] for _ in range(ny)], dtype=np.float32)
    obstacles = np.zeros((ny, nx), dtype=np.bool_)
    for x, y in block:
            obstacles[x,y] = 1
    av_vels = [0.0] * param["maxIters"]

    logger.debug(param)

    init_toc = time.time_ns()
    comp_tic = time.time_ns()

    for tt in range(param["maxIters"]):
        # logger.info("==Start step==")
        timestep(cells, obstacles, nx, ny, density, accel, omega)
        # av_vels[tt] = av_velocity(cells, obstacles, nx, ny)
        logger.info(f"==timestep: {tt}==")
        # logger.info(f"ev velocity: {av_vels[tt]}")
        # logger.info(f"tot density: {total_density(cells, nx, ny)}")

    comp_toc = time.time_ns()
    col_tic = time.time_ns()

    # TODO collate data here

    col_toc = time.time_ns()
    tot_toc = time.time_ns()

    print(f"Reynolds number:\t\t{calc_reynolds(cells, obstacles, nx, ny, omega, reynolds_dim)}")
    print(f"Elapsed Init time:\t\t\t{init_toc - init_tic}")
    print(f"Elapsed Compute time:\t\t\t{comp_toc - comp_tic}")
    print(f"Elapsed Collate time:\t\t\t{col_toc - col_tic}")
    print(f"Elapsed Total time:\t\t\t{tot_toc - tot_tic}")

    write_values(cells, obstacles, av_vels, nx, ny, density, param["maxIters"])


def calc_reynolds(cells, obstacles, nx, ny, omega, reynolds_dim):
    viscosity = 1. / 6. * (2. / omega - 1.)
    return av_velocity(cells, obstacles, nx, ny) * reynolds_dim / viscosity


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
                    local_density = 0.
                    for kk in range(NSPEEDS):
                        local_density += cells[ii][jj][kk]
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
            f.write(f"{ii} {jj} {u_x} {u_y} {u} {pressure} {obstacles[ii, jj]}")

    with open(AVVELSFILE, "w") as f:
        for ii in range(maxIters):
            f.write(f"{ii} {av_vels[ii]}")


def av_velocity(cells, obstacles, nx, ny):
    """TODO"""
    tot_cells = 0
    tot_u = 0.
    for jj in range(ny):
        for ii in range(nx):
            if not obstacles[ii][jj]:
                local_density = 0.
                for kk in range(NSPEEDS):
                    local_density += cells[ii][jj][kk]
                # x-component of velocity
                u_x = (
                              cells[ii][jj][1] + cells[ii][jj][5] +
                              cells[ii][jj][8] - (
                                      cells[ii][jj][3] + cells[ii][jj][6] +
                                      cells[ii][jj][7]
                              )
                      ) / local_density
                # compute y velocity component
                u_y = (
                              cells[ii][jj][2] + cells[ii][jj][5] +
                              cells[ii][jj][6] - (
                                      cells[ii][jj][4] + cells[ii][jj][7] +
                                      cells[ii][jj][8]
                              )
                      ) / local_density
                # accumulate the norm of x- and y- velocity components
                tot_u += sqrt((u_x * u_x) + (u_y * u_y))
                # increase counter of inspected cells
                tot_cells += 1
    return tot_u / tot_cells

def total_density(cells, nx, ny):
    total = 0.
    for jj in range(ny):
        for ii in range(nx):
            for kk in range(NSPEEDS):
                total += cells[ii][jj][kk]
    return total

@njit(parallel=True)
def timestep(cells, obstacles, nx, ny, density, accel, omega):
    """One step elapse"""
    c_sq = 1. / 3.  # square of speed of sound
    w0 = 4. / 9.  # weighting factor
    w1 = 1. / 9.  # weighting factor
    w2 = 1. / 36.  # weighting factor
    w3 = density * accel / 9.
    w4 = density * accel / 36.

    tmp = cells.copy()
    c = np.zeros(NSPEEDS, dtype=np.float32)

    u = np.zeros((NSPEEDS), dtype=np.float32)
    d_equ = np.zeros((NSPEEDS), dtype=np.float32)

    for jj in range(ny):
        for ii in range(nx):
            # accelerate and propagate
            y_n = (jj + 1) % ny
            x_e = (ii + 1) % nx
            y_s = (ny - 1) if jj == 0 else jj - 1
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

if __name__ == "__main__":
    main()
