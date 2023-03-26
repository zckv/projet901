from argparse import ArgumentParser
from math import sqrt

import coloredlogs
import logging
import time
import os
import sys

from sys import argv

NSPEEDS = 9
FINALSTATEFILE = "final_state.dat"
INITIALSTATEFILE = "initial_state.dat"
AVVELSFILE = "av_vels.dat"

logging.basicConfig(encoding='utf-8', level=logging.DEBUG)
logger = logging.getLogger(os.path.basename(sys.argv[0]))
logger.setLevel(logging.DEBUG)
coloredlogs.DEFAULT_FIELD_STYLES["levelname"]["color"] = "cyan"
coloredlogs.install(logger=logger, level=logging.DEBUG)


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
            block.add((x, y))
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
    cells = [[[w0, w1, w1, w1, w1, w2, w2, w2, w2] for _ in range(nx)] for _ in range(ny)]
    tmp_cells = [[[0.0] * NSPEEDS for _ in range(nx)] for _ in range(ny)]
    obstacles = [[1 if (x, y) in block else 0 for x in range(nx)] for y in range(ny)]
    av_vels = [0.0] * param["maxIters"]

    logger.debug(param)

    init_toc = time.time_ns()
    comp_tic = time.time_ns()

    for tt in range(param["maxIters"]):
        logger.debug("==Start step==")
        timestep(cells, tmp_cells, obstacles, nx, ny, density, accel, omega)
        av_vels[tt] = av_velocity(cells, obstacles, nx, ny)
        logger.debug(f"==timestep: {tt}==")
        logger.debug(f"ev velocity: {av_vels[tt]}")
        logger.debug(f"tot density: {total_density(cells, nx, ny)}")

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
                if (ii, jj) in obstacles:
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
            f.write(f"{ii} {jj} {u_x} {u_y} {u} {pressure} {1 if (ii, jj) in obstacles else 0}")

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


def timestep(cells, tmp_cells, obstacles, nx, ny, density, accel, omega):
    """One step elapse"""
    accelerate_flow(cells, obstacles, density, accel, nx, ny)
    propagate(cells, tmp_cells, nx, ny)
    rebound(cells, tmp_cells, obstacles, nx, ny)
    collision(cells, tmp_cells, obstacles, nx, ny, omega)


def accelerate_flow(cells, obstacles, density, accel, nx, ny):
    # compute weighting factors
    w1 = density * accel / 9.
    w2 = density * accel / 36.

    # modify the 2nd row of the grid
    jj = ny - 2
    for ii in range(nx):
        # if the cell is not occupied, and we don't send a negative density
        if (ii, jj) not in obstacles and \
                (cells[ii][jj][3] - w1) > 0 and \
                (cells[ii][jj][6] - w2) > 0 and \
                (cells[ii][jj][7] - w2) > 0:
            # increase 'east-side' densities
            cells[ii][jj][1] += w1
            cells[ii][jj][5] += w2
            cells[ii][jj][8] += w2
            # decrease 'west-side' densities
            cells[ii][jj][3] -= w1
            cells[ii][jj][6] -= w2
            cells[ii][jj][7] -= w2


def propagate(cells, tmp_cells, nx, ny):
    """TODO Optimize and parallelize"""
    for jj in range(ny):
        for ii in range(nx):
            y_n = (jj + 1) % ny
            x_e = (ii + 1) % nx
            y_s = (ny - 1) if jj == 0 else jj - 1
            x_w = (nx - 1) if ii == 0 else ii - 1

            # propagate densities from neighbouring cells, following
            # appropriate directions of  travel and writing into
            # scratch space grid

            tmp_cells[ii][jj][0] = cells[ii][jj][0]
            tmp_cells[ii][jj][1] = cells[x_e][jj][1]
            tmp_cells[ii][jj][2] = cells[ii][y_n][2]
            tmp_cells[ii][jj][3] = cells[x_w][jj][3]
            tmp_cells[ii][jj][4] = cells[ii][y_s][4]
            tmp_cells[ii][jj][5] = cells[x_e][y_n][5]
            tmp_cells[ii][jj][6] = cells[x_w][y_n][6]
            tmp_cells[ii][jj][7] = cells[x_w][y_s][7]
            tmp_cells[ii][jj][8] = cells[x_e][y_s][8]


def rebound(cells, tmp_cells, obstacles, nx, ny):
    for jj in range(ny):
        for ii in range(nx):
            if (ii, jj) in obstacles:
                # called after propagate, so taking values from scratch space
                cells[ii][jj][1] = tmp_cells[ii][jj][3]
                cells[ii][jj][2] = tmp_cells[ii][jj][4]
                cells[ii][jj][3] = tmp_cells[ii][jj][1]
                cells[ii][jj][4] = tmp_cells[ii][jj][2]
                cells[ii][jj][5] = tmp_cells[ii][jj][7]
                cells[ii][jj][6] = tmp_cells[ii][jj][8]
                cells[ii][jj][7] = tmp_cells[ii][jj][5]
                cells[ii][jj][8] = tmp_cells[ii][jj][6]


def collision(cells, tmp_cells, obstacles, nx, ny, omega):
    """TODO Optimize and parallelize
    loop over the cells in the grid
    NB the collision step is called after
    the propagate step and so values of interest
    are in the scratch-space grid"""

    c_sq = 1. / 3.  # square of speed of sound
    w0 = 4. / 9.  # weighting factor
    w1 = 1. / 9.  # weighting factor
    w2 = 1. / 36.  # weighting factor

    # pragma omp parallel for
    for jj in range(ny):
        # pragma omp parallel for
        for ii in range(nx):
            if not obstacles[ii][jj]:
                # compute local density total
                local_density = 0.

                # pragma omp parallel for
                for kk in range(NSPEEDS):
                    local_density += tmp_cells[ii][jj][kk]

                # compute x velocity component
                u_x = (
                              tmp_cells[ii][jj][1]
                              + tmp_cells[ii][jj][5]
                              + tmp_cells[ii][jj][8]
                              - (
                                      tmp_cells[ii][jj][3]
                                      + tmp_cells[ii][jj][6]
                                      + tmp_cells[ii][jj][7]
                              )
                      ) / local_density

                # compute y velocity component
                u_y = (
                              tmp_cells[ii][jj][2]
                              + tmp_cells[ii][jj][5]
                              + tmp_cells[ii][jj][6]
                              - (
                                      tmp_cells[ii][jj][4]
                                      + tmp_cells[ii][jj][7]
                                      + tmp_cells[ii][jj][8]
                              )
                      ) / local_density

                # velocity squared
                u_sq = u_x * u_x + u_y * u_y

                # directional velocity components
                u = [0] * NSPEEDS
                u[1] = u_x
                u[2] = u_y
                u[3] = -u_x
                u[4] = -u_y
                u[5] = u_x + u_y
                u[6] = -u_x + u_y
                u[7] = -u_x - u_y
                u[8] = u_x - u_y

                # equilibrium densities
                d_equ = [0] * NSPEEDS
                # zero velocity density: weight w0
                d_equ[0] = w0 * local_density * (1. - u_sq / (2. * c_sq))

                # axis speeds: weight w1 */
                d_equ[1] = w1 * local_density * (
                        1 + u[1] / c_sq + (u[1] * u[1]) / (2 * c_sq * c_sq) - u_sq / (2 * c_sq))
                d_equ[2] = w1 * local_density * (
                        1 + u[2] / c_sq + (u[2] * u[2]) / (2 * c_sq * c_sq) - u_sq / (2 * c_sq))
                d_equ[3] = w1 * local_density * (
                        1 + u[3] / c_sq + (u[3] * u[3]) / (2 * c_sq * c_sq) - u_sq / (2 * c_sq))
                d_equ[4] = w1 * local_density * (
                        1 + u[4] / c_sq + (u[4] * u[4]) / (2 * c_sq * c_sq) - u_sq / (2 * c_sq))

                # diagonal speeds: weight w2 */
                d_equ[5] = w2 * local_density * (
                        1 + u[5] / c_sq + (u[5] * u[5]) / (2 * c_sq * c_sq) - u_sq / (2 * c_sq))
                d_equ[6] = w2 * local_density * (
                        1 + u[6] / c_sq + (u[6] * u[6]) / (2 * c_sq * c_sq) - u_sq / (2 * c_sq))
                d_equ[7] = w2 * local_density * (
                        1 + u[7] / c_sq + (u[7] * u[7]) / (2 * c_sq * c_sq) - u_sq / (2 * c_sq))
                d_equ[8] = w2 * local_density * (
                        1 + u[8] / c_sq + (u[8] * u[8]) / (2 * c_sq * c_sq) - u_sq / (2 * c_sq))

                # relaxation step
                for kk in range(NSPEEDS):
                    cells[ii][jj][kk] = tmp_cells[ii][jj][kk] + omega * (
                            d_equ[kk] - tmp_cells[ii][jj][kk])


if __name__ == "__main__":
    main()
