

def timestep(cells, tmp, obstacles, nx, ny, density, accel, omega):
    c_sq = 1. / 3.  # square of speed of sound
    w0 = 4. / 9.  # weighting factor
    w1 = 1. / 9.  # weighting factor
    w2 = 1. / 36.  # weighting factor
    w3 = density * accel / 9.
    w4 = density * accel / 36.

    for ii in range(nx):
        for jj in range(ny):
            # soit
            # si pas d'obstacle
            if not obstacles[ii][jj]:
                z =  np.array(
                    cells[ii][jj][0],
                    cells[x_e][jj][1],
                    cells[ii][y_n][2],
                    cells[x_w][jj][3],
                    cells[ii][y_s][4],
                    cells[x_e][y_n][5],
                    cells[x_w][y_n][6],
                    cells[x_w][y_s][7],
                    cells[x_e][y_s][8],
                )
                s_z = sum(z)

                u_x = (z[1] + z[5] + z[8] - z[3] - z[6] - z[7]) / s_z
                u_y = (z[2] + z[5] + z[6] - z[4] - z[7] - z[8]) / s_z
                u_sq = u_x * u_x + u_y * u_y

                cells[ii][jj][0] = z[0] + omega * (
                    (
                        w0 * sum(z) * (1. - u_sq / (2. * c_sq))
                    ) - z[0]
                )

                cells[ii][jj][1] = z[1] + omega * (
                    (
                        w1 * s_z * (
                            1 + u_x/c_sq +
                            u_x**2/(2*c_sq*c_sq) -
                            u_sq/(2*c_sq)
                        )
                    ) - z[1]
                )

                cells[ii][jj][2] = z[2] + omega * (
                    (
                        w1 * sum(z) * (1. - u_sq / (2. * c_sq))
                    ) - z[2]
                )


