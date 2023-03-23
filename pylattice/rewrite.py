from argparse import ArgumentParser
from math import sqrt
import logging

NSPEEDS = 9
FINALSTATEFILE = "final_state.dat"
INITIALSTATEFILE = "initial_state.dat"
AVVELSFILE = "av_vels.dat"

logger = logging.Logger(name="MainLogger")


class ProjectParser(ArgumentParser):
    """Parse args from command line"""
    # TODO


def read_input_file(file_path):
    d = {}
    with open(file_path, 'r') as f:
        d["nx"] = f.readLine()
        d["ny"] = f.readLine()
        d["maxIters"] = f.readLine()
        d["reynolds_dim"] = f.readLine()
        d["density"] = f.readLine()
        d["accel"] = f.readLine()
        d["omega"] = f.readLine()
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
    param = read_input_file("")
    block = read_obstacle_file("")

    w0 = param["density"] * 4 / 9
    w1 = param["density"] / 9
    w2 = param["density"] / 36

    nx = param["nx"]
    ny = param["ny"]
    density = param["density"]
    accel = param["accel"]
    omega = param["omega"]

    # main_grid
    cells = [[[w0, w1, w1, w1, w1, w2, w2, w2, w2] for _ in range(nx)] for _ in range(ny)]
    tmp_cells = [[[0.0 * NSPEEDS] for _ in range(nx)] for _ in range(ny)]
    obstacles = [[1 if (x, y) in block else 0 for x in range(nx)] for y in range(ny)]
    av_vels = [0.0] * param["maxIters"]

    logger.debug(param)

    for tt in range(param["maxIters"]):
        timestep(cells, tmp_cells, obstacles, nx, ny, density, accel, omega)
        av_vels[tt] = av_velocity(cells, obstacles, nx, ny)
        logger.debug(f"==timestep: {tt}==")
        logger.debug(f"ev velocity: {av_vels[tt]}==")
        logger.debug(f"tot density: {total_density(cells, nx, ny)}==")


def av_velocity(cells, obstacles, nx, ny):
    """TODO"""
    tot_cells = 0
    tot_u = 0.

    for jj in range(ny):
        for ii in range(nx):
            if not obstacles[ii + jj * nx]:
                local_density = 0.
                for kk in range(NSPEEDS):
                    local_density += cells[ii + jj * nx].speeds[kk]

                # x-component of velocity
                u_x = (
                              cells[ii + jj * nx].speeds[1] + cells[ii + jj * nx].speeds[5] +
                              cells[ii + jj * nx].speeds[8] - (
                                      cells[ii + jj * nx].speeds[3] + cells[ii + jj * nx].speeds[6] +
                                      cells[ii + jj * nx].speeds[7]
                              )
                      ) / local_density
                # compute y velocity component

                u_y = (
                              cells[ii + jj * nx].speeds[2] + cells[ii + jj * nx].speeds[5] +
                              cells[ii + jj * nx].speeds[6] - (
                                      cells[ii + jj * nx].speeds[4] + cells[ii + jj * nx].speeds[7] +
                                      cells[ii + jj * nx].speeds[8]
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
                total += cells[ii + jj * nx].speeds[kk]


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
        if (not obstacles[ii + jj * nx]) and \
                (cells[ii + jj * nx][3] - w1) > 0 and \
                (cells[ii + jj * nx][6] - w2) > 0 and \
                (cells[ii + jj * nx][7] - w2) > 0:
            # increase 'east-side' densities
            cells[ii + jj * nx][1] += w1
            cells[ii + jj * nx][5] += w2
            cells[ii + jj * nx][8] += w2
            # decrease 'west-side' densities
            cells[ii + jj * nx][3] -= w1
            cells[ii + jj * nx][6] -= w2
            cells[ii + jj * nx][7] -= w2


def propagate(cells, tmp_cells, nx, ny):
    """TODO Optimize and parallelize"""
    for jj in range(ny):
        for ii in range(nx):
            y_n = (jj + 1) % ny
            x_e = (ii + 1) % nx
            y_s = (jj + ny - 1) if (jj == 0) else (jj - 1)
            x_w = (ii + nx - 1) if (ii == 0) else (ii - 1)
            # propagate densities from neighbouring cells, following
            # appropriate directions of  travel and writing into
            # scratch space grid

            tmp_cells[ii + jj * nx][0] = cells[ii + jj * nx][0]
            tmp_cells[ii + jj * nx][1] = cells[x_w + jj * nx][1]
            tmp_cells[ii + jj * nx][2] = cells[ii + y_s * nx][2]
            tmp_cells[ii + jj * nx][3] = cells[x_e + jj * nx][3]
            tmp_cells[ii + jj * nx][4] = cells[ii + y_n * nx][4]
            tmp_cells[ii + jj * nx][5] = cells[x_w + y_s * nx][5]
            tmp_cells[ii + jj * nx][6] = cells[x_e + y_s * nx][6]
            tmp_cells[ii + jj * nx][7] = cells[x_e + y_n * nx][7]
            tmp_cells[ii + jj * nx][8] = cells[x_w + y_n * nx][8]


def rebound(cells, tmp_cells, obstacles, nx, ny):
    for jj in range(ny):
        for ii in range(nx):
            if obstacles[jj*nx + ii]:
                # called after propagate, so taking values from scratch space
                cells[ii + jj * nx][1] = tmp_cells[ii + jj * nx][3]
                cells[ii + jj * nx][2] = tmp_cells[ii + jj * nx][4]
                cells[ii + jj * nx][3] = tmp_cells[ii + jj * nx][1]
                cells[ii + jj * nx][4] = tmp_cells[ii + jj * nx][2]
                cells[ii + jj * nx][5] = tmp_cells[ii + jj * nx][7]
                cells[ii + jj * nx][6] = tmp_cells[ii + jj * nx][8]
                cells[ii + jj * nx][7] = tmp_cells[ii + jj * nx][5]
                cells[ii + jj * nx][8] = tmp_cells[ii + jj * nx][6]


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
            if not obstacles[ii + jj * nx]:
                # compute local density total
                local_density = 0.

                # pragma omp parallel for
                for kk in range(NSPEEDS):
                    local_density += tmp_cells[ii + jj * nx][kk]

                # compute x velocity component
                u_x = (
                              tmp_cells[ii + jj * nx][1]
                              + tmp_cells[ii + jj * nx][5]
                              + tmp_cells[ii + jj * nx][8]
                              - (
                                      tmp_cells[ii + jj * nx][3]
                                      + tmp_cells[ii + jj * nx][6]
                                      + tmp_cells[ii + jj * nx][7]
                              )
                      ) / local_density

                # compute y velocity component
                u_y = (
                              tmp_cells[ii + jj * nx][2]
                              + tmp_cells[ii + jj * nx][5]
                              + tmp_cells[ii + jj * nx][6]
                              - (
                                      tmp_cells[ii + jj * nx][4]
                                      + tmp_cells[ii + jj * nx][7]
                                      + tmp_cells[ii + jj * nx][8]
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
                    cells[ii + jj * nx][kk] = tmp_cells[ii + jj * nx][kk] + omega * (
                            d_equ[kk] - tmp_cells[ii + jj * nx][kk])


if __name__ == "__main__":
    pass