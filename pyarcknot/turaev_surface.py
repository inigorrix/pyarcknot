import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


# Get crossing surface points
def crossing_surfpoints(base_size=20, height=10, curve=5, res_c=4, res_z=6):
    b, h, c = base_size, height, curve
    nc, nz = res_c, round(res_z * .5) + 1
    ta, curve = np.linspace(0, 1, nc), np.empty([2, nc])
    zh = np.linspace(h / 2, 0, nz)
    xtot, ytot = np.empty([2 * nz, nc]), np.empty([2 * nz, nc])
    ztot = np.empty([2 * nz, nc])

    P0 = np.array([0, b / 2])
    for j, z in enumerate(zh):
        c1 = (z * c) ** .5
        c2 = c1 * .5
        P1, P2 = np.array([0, c1]), np.array([c2, c2])
        for i in range(2):
            curve[i] = (1 - ta) ** 2 * P0[i] + 2 * (1 - ta) * ta * P1[i] + ta ** 2 * P2[i]

        k = -(j + 1)
        xtot[j], ytot[j], ztot[j] = curve[0], curve[1], np.linspace(z, z, nc)
        xtot[k], ytot[k], ztot[k] = -curve[0], curve[1], -np.linspace(z, z, nc)

    return xtot, ytot, ztot


# Draw crossings surfaces
def draw_crossing(ax, smoothing=4, coords=(10, 10),
                  colormap=cm.viridis, surfpoints=crossing_surfpoints()):
    xq, yq, zq = surfpoints
    if smoothing == 5:
        xq = -xq
    x_c, y_c = coords[0], coords[1]
    ax.plot_surface(xq + x_c, yq + y_c, zq, cmap=colormap)
    ax.plot_surface(yq + x_c, xq + y_c, zq, cmap=colormap)
    ax.plot_surface(-yq + x_c, -xq + y_c, zq, cmap=colormap)
    ax.plot_surface(-xq + x_c, -yq + y_c, zq, cmap=colormap)


# Draw lines surfaces
def draw_line(ax, coords0, coords1, height=10, res_z=6, colormap=cm.viridis):
    h, nz = height * .5, round(res_z * .5) * 2 + 1
    x_coords = np.array([coords0[0], coords1[0]])
    y_coords = np.array([coords0[1], coords1[1]])
    z_coords = np.zeros([nz, 2])
    z_coords[:, 0], z_coords[:, 1] = np.linspace(h, -h, nz), np.linspace(h, -h, nz)

    ax.plot_surface(x_coords, y_coords, z_coords, cmap=colormap)


def draw_element(element, curve=5, res_c=4, res_z=6,
                 coords0=(0, 0), coords1=(10, 10),
                 axis_on=True, view=(-30, 30)):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax._axis3don = axis_on
    if element == "crossing":
        draw_crossing(ax, surfpoints=crossing_surfpoints(curve=curve,
                                                         res_c=res_c,
                                                         res_z=res_z))
    elif element == "line":
        draw_line(ax, coords0, coords1, res_z=res_z)
    
    set_3daxes_equal(ax)
    ax.view_init(azim=view[0], elev=view[1])
    plt.show()


def set_3daxes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


# Draw Turaev surface of a diagram
def turaev_surf(xabo, figsize=None, size=20, height=10, res_z=6,
                view=(-45, 30), colormap=cm.viridis):
    if figsize is not None:
        figsize = (figsize, figsize)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection='3d')
    ax._axis3don = False

    s = np.shape(xabo)[0]

    ## VERTICAL AND HORIZONTAL LINES
    for i in range(s):
        row, col = xabo[i, :], xabo[:, i]

        hverts, n_hlines = np.where(row > 0)[0], np.count_nonzero(row > 0) - 1
        for j in range(n_hlines):
            c0, c1 = hverts[j:j + 2]
            if row[c0] > 2: c0 += .5
            if row[c1] > 2: c1 -= .5
            p0, p1 = np.array([c0, s - (i + 1)]) * size, np.array([c1, s - (i + 1)]) * size
            draw_line(ax, p0, p1, height=height, res_z=res_z,
                      colormap=colormap)

        vverts, n_vlines = np.where(col > 0)[0], np.count_nonzero(col > 0) - 1
        for j in range(n_vlines):
            r0, r1 = vverts[j:j + 2]
            if j != 0: r0 += .5
            if j != n_vlines - 1: r1 -= .5
            p0, p1 = np.array([i, s - (r0 + 1)]) * size, np.array([i, s - (r1 + 1)]) * size
            draw_line(ax, p0, p1, height=height, res_z=res_z,
                      colormap=colormap)

    ## CROSSINGS
    cr_surfpoints = crossing_surfpoints(base_size=size, height=height,
                                        res_z=res_z)
    c_pos, n_c = np.where(xabo > 2), np.count_nonzero(xabo > 2)
    for i in range(n_c):
        row, col = c_pos[0][i], c_pos[1][i]
        ab = xabo[row, col]
        cr_coords = np.array([col, s - (row + 1)]) * size
        draw_crossing(ax, ab, cr_coords,
                      surfpoints=cr_surfpoints, colormap=colormap)

    ax.view_init(azim=view[0], elev=view[1])
    set_3daxes_equal(ax)
    plt.show()

