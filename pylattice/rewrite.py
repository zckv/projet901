from argparse import ArgumentParser
import logging
NSPEEDS = 9

logger = logging.Logger()


class ProjectParser(ArgumentParser):
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
            block.add((x,y))
    return block

def main():
    # TODO initalization
    param = read_input_file("")
    block = read_obstacle_file("")

    w0 = param["density"] * 4 / 9
    w1 = param["density"] / 9
    w2 = param["density"] / 36

    # main_grid
    cells = [[[w0, w1, w1, w1, w1, w2, w2, w2, w2] for _ in range(nx)] for _ in range(ny)]
    tmp_cells = [[[0.0 * NSPEEDS] for _ in range(nx)] for _ in range(ny)]
    obstacles = [[1 if (x,y) in block else 0 for x in range(nx)] for y in range(ny)]
    av_vels = [0.0] * param["maxIters"]

    logger.debug(param)

    for tt in range(param["maxIters"]):
        timestep()
        av_vels[tt] = av_velocity()
        logger.debug(av_vels)

def av_velocity():
    pass

def timestep():
    pass
