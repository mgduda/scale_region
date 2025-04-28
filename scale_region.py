#!/usr/bin/env python3

def stereo_project(lam, phi, lam_0, phi_1):
    from math import sin, cos
    R = 1.0
    k = 2.0 * R / (1.0 + sin(phi_1) * sin(phi)
                   + cos(phi_1) * cos(phi) * cos(lam - lam_0))
    x = k * cos(phi) * sin(lam - lam_0)
    y = k * (cos(phi_1) * sin(phi)
             - sin(phi_1) * cos(phi) * cos(lam - lam_0))
    return (x, y)


def stereo_inverse(x, y, lam_0, phi_1):
    from math import sqrt, asin, sin, cos, atan2
    R = 1.0
    rho = sqrt(x * x + y * y)
    c = 2.0 * atan2(rho, 2.0 * R)
    phi = (asin(cos(c) * sin(phi_1)
           + y * sin(c) * cos(phi_1) / rho))
    lam = lam_0 + atan2(x * sin(c),
                             rho * cos(phi_1) * cos(c)
                             - y * sin(phi_1) * sin(c))
    return (lam, phi)


def geo_to_cart(lam, phi):
    from math import cos, sin
    R = 1.0
    z = R * sin(phi)
    x = R * cos(lam) * cos(phi)
    y = R * sin(lam) * cos(phi)
    return (x, y, z)


def cart_to_geo(x, y, z):
    from math import atan2, asin
    lam = atan2(y, x)
    phi = asin(z)
    return (lam, phi)


def scale_var(meshfile, varname, scalefac):
    tmp = meshfile.variables[varname][:]
    meshfile.variables[varname][:] = tmp / scalefac


def scale_location(meshfile, loc, scalefac, cenlon, cenlat):
    lonLoc = meshfile.variables['lon'+loc][:]
    latLoc = meshfile.variables['lat'+loc][:]
    xLoc = meshfile.variables['x'+loc][:]
    yLoc = meshfile.variables['y'+loc][:]
    zLoc = meshfile.variables['z'+loc][:]

    for i in range(len(lonLoc)):
        (x,y) = stereo_project(lonLoc[i], latLoc[i], cenlon, cenlat)
        x = x / scalefac
        y = y / scalefac
        (lam, phi) = stereo_inverse(x, y, cenlon, cenlat)
        lonLoc[i] = lam
        latLoc[i] = phi
        (x, y, z) = geo_to_cart(lam, phi)
        xLoc[i] = x
        yLoc[i] = y
        zLoc[i] = z

    meshfile.variables['lon'+loc][:] = lonLoc
    meshfile.variables['lat'+loc][:] = latLoc
    meshfile.variables['x'+loc][:] = xLoc
    meshfile.variables['y'+loc][:] = yLoc
    meshfile.variables['z'+loc][:] = zLoc


def plane_distance(a, b):
    from math import sqrt

    d = b - a

    return sqrt(np.sum(d * d))


def sphere_distance(a, b):
    from math import asin

    c = plane_distance(a, b)

    return 2.0 * asin(c / 2.0)


def triangle_area(A, B, C):
    from math import sqrt, tan, atan

    a = sphere_distance(B, C)
    b = sphere_distance(A, C)
    c = sphere_distance(A, B)

    # Compute semi-perimeter
    s = 0.5 * (a + b + c)

    tanqe = sqrt(max(0.0,
                     tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c))))

    return 4.0 * atan(tanqe)


def sign_for_edge(cellID, edgeID, cellsOnEdge):
    if cellID == cellsOnEdge[edgeID, 0]:
        return 1.0
    else:
        return -1.0


def rotate_kiteAreasOnCell(kiteAreasOnCell, edgeID, edgesOnCell, nEdgesOnCell):
    rotated_kiteAreasOnCell = np.zeros(nEdgesOnCell)

    #
    # Find starting edge in the edgesOnCell list
    # NB: verticesOnCell (and therefore, kiteAreasOnCell), "lead" edgesOnCell,
    #     so we actually need to start with the vertex whose index is one more
    #     than the index of the edgeID
    #
    for i in range(nEdgesOnCell):
        if edgesOnCell[i] == edgeID:
            istart = (i+1) % nEdgesOnCell
            break

    #
    # Copy kiteAreasOnCell
    #
    j = 0
    for i in range(istart, nEdgesOnCell):
        rotated_kiteAreasOnCell[j] = kiteAreasOnCell[i]
        j = j + 1

    for i in range(istart):
        rotated_kiteAreasOnCell[j] = kiteAreasOnCell[i]
        j = j + 1

    return rotated_kiteAreasOnCell


def find_neighbor(cellsOnCell, maxEdges, neighbor):
    for i in range(maxEdges):
        if cellsOnCell[i] == neighbor:
            return i

    return -1


if __name__ == '__main__':
    import argparse
    import os
    import sys
    import shutil
    import numpy as np
    from netCDF4 import Dataset

    parser = argparse.ArgumentParser()
    parser.add_argument('orig_mesh',
                        help='the original regional mesh netCDF file')
    parser.add_argument('scaled_mesh',
                        help='the scaled regional mesh netCDF file to be created')
    parser.add_argument('scale_fac', type=float,
                        help='the scale factor. Values >1.0 increase the mesh resolution.')
    parser.add_argument('tan_lat', type=float,
                        help='the latitude of the tangent point, in degrees')
    parser.add_argument('tan_lon', type=float,
                        help='the longitude of the tangent point, in degrees')
    args = parser.parse_args()

    #
    # Begin by creating a copy of the original regional mesh
    #
    if not os.path.isfile(args.orig_mesh):
        print('')
        print('Error: Could not find original regional mesh file '
              + args.orig_mesh)
        print('')
        sys.exit(1)

    try:
        print('Copying ' + args.orig_mesh + ' to ' + args.scaled_mesh)
        shutil.copy(args.orig_mesh, args.scaled_mesh)
    except:
        print('')
        print('Error: Could not copy ' + args.orig_mesh + ' to '
              + args.scaled_mesh)
        print('')
        sys.exit(1)


    #
    # Scale mesh fields in the new netCDF file
    #
    try:
        mesh = Dataset(args.scaled_mesh,'r+')
    except:
        print('')
        print('Error: Could not open ' + args.scaled_mesh + ' as a netCDF file')
        print('')
        sys.exit(1)

    nCells = mesh.dimensions['nCells'].size
    nVertices = mesh.dimensions['nVertices'].size
    nEdges = mesh.dimensions['nEdges'].size
    vertexDegree = mesh.dimensions['vertexDegree'].size
    maxEdges = mesh.dimensions['maxEdges'].size
    maxEdges2 = mesh.dimensions['maxEdges2'].size

    scalefac = args.scale_fac
    cenlat = args.tan_lat
    cenlon = args.tan_lon

    print('Recomputing cell coordinates...')
    scale_location(mesh, 'Cell', scalefac, cenlon, cenlat)

    print('Recomputing vertex coordinates...')
    scale_location(mesh, 'Vertex', scalefac, cenlon, cenlat)

    print('Recomputing edge coordinates...')
    scale_location(mesh, 'Edge', scalefac, cenlon, cenlat)

    scale_var(mesh, 'nominalMinDc', scalefac)


    #
    # Recompute distances and areas after scaling
    #
    xCell = mesh.variables['xCell'][:]
    yCell = mesh.variables['yCell'][:]
    zCell = mesh.variables['zCell'][:]
    xVertex = mesh.variables['xVertex'][:]
    yVertex = mesh.variables['yVertex'][:]
    zVertex = mesh.variables['zVertex'][:]
    xEdge = mesh.variables['xEdge'][:]
    yEdge = mesh.variables['yEdge'][:]
    zEdge = mesh.variables['zEdge'][:]

    nEdgesOnCell = mesh.variables['nEdgesOnCell'][:]
    nEdgesOnEdge = mesh.variables['nEdgesOnEdge'][:]

    cellsOnEdge = mesh.variables['cellsOnEdge'][:] - 1
    verticesOnEdge = mesh.variables['verticesOnEdge'][:] - 1
    edgesOnEdge = mesh.variables['edgesOnEdge'][:] - 1
    cellsOnVertex = mesh.variables['cellsOnVertex'][:] - 1
    edgesOnVertex = mesh.variables['edgesOnVertex'][:] - 1
    verticesOnCell = mesh.variables['verticesOnCell'][:] - 1
    edgesOnCell = mesh.variables['edgesOnCell'][:] - 1

    dcEdge = mesh.variables['dcEdge'][:]
    dvEdge = mesh.variables['dvEdge'][:]

    areaTriangle = mesh.variables['areaTriangle'][:]
    areaCell = mesh.variables['areaCell'][:]
    kiteAreasOnVertex = mesh.variables['kiteAreasOnVertex'][:]

    weightsOnEdge = mesh.variables['weightsOnEdge'][:]

    #
    # dcEdge
    #
    print('Recomputing dcEdge...')
    dcEdge_new = np.zeros((nEdges))
    for i in range(nEdges):
        coe1 = cellsOnEdge[i, 0]
        coe2 = cellsOnEdge[i, 1]
        if coe1 >= 0 and coe2 >= 0:
            cell1 = np.asarray([xCell[coe1], yCell[coe1], zCell[coe1]])
            cell2 = np.asarray([xCell[coe2], yCell[coe2], zCell[coe2]])
            dcEdge_new[i] = sphere_distance(cell1, cell2)
        else:
            dcEdge_new[i] = dcEdge[i] / scalefac

    mesh.variables['dcEdge'][:] = dcEdge_new
    dcEdge = dcEdge_new

    #
    # dvEdge
    #
    print('Recomputing dvEdge...')
    dvEdge_new = np.zeros((nEdges))
    for i in range(nEdges):
        voe1 = verticesOnEdge[i, 0]
        voe2 = verticesOnEdge[i, 1]
        vtx1 = np.asarray([xVertex[voe1], yVertex[voe1], zVertex[voe1]])
        vtx2 = np.asarray([xVertex[voe2], yVertex[voe2], zVertex[voe2]])
        dvEdge_new[i] = sphere_distance(vtx1, vtx2)

    mesh.variables['dvEdge'][:] = dvEdge_new
    dvEdge = dvEdge_new

    #
    # areaTriangle
    #
    print('Recomputing areaTriangle...')
    areaTriangle_new = np.zeros((nVertices))
    for i in range(nVertices):
        cov1 = cellsOnVertex[i, 0]
        cov2 = cellsOnVertex[i, 1]
        cov3 = cellsOnVertex[i, 2]
        if cov1 >= 0 and cov2 >= 0 and cov3 >= 0:
            cell1 = np.asarray([xCell[cov1], yCell[cov1], zCell[cov1]])
            cell2 = np.asarray([xCell[cov2], yCell[cov2], zCell[cov2]])
            cell3 = np.asarray([xCell[cov3], yCell[cov3], zCell[cov3]])
            areaTriangle_new[i] = triangle_area(cell1, cell2, cell3)
        else:
            areaTriangle_new[i] = areaTriangle[i] / scalefac**2

    mesh.variables['areaTriangle'][:] = areaTriangle_new
    areaTriangle = areaTriangle_new

    #
    # areaCell
    #
    print('Recomputing areaCell...')
    areaCell_new = np.zeros((nCells))
    for i in range(nCells):
        cell = np.asarray([xCell[i], yCell[i], zCell[i]])
        for j in range(nEdgesOnCell[i]):
            voc1 = verticesOnCell[i, j]
            voc2 = verticesOnCell[i, (j+1)%nEdgesOnCell[i]]
            vtx1 = np.asarray([xVertex[voc1], yVertex[voc1], zVertex[voc1]])
            vtx2 = np.asarray([xVertex[voc2], yVertex[voc2], zVertex[voc2]])
            areaCell_new[i] = areaCell_new[i] + triangle_area(cell, vtx1, vtx2)

    mesh.variables['areaCell'][:] = areaCell_new
    areaCell = areaCell_new

    #
    # kiteAreasOnVertex
    #
    print('Recomputing kiteAreasOnVertex...')
    kiteAreasOnVertex_new = np.zeros((nVertices,vertexDegree))
    for i in range(nVertices):
        vtx = np.asarray([xVertex[i], yVertex[i], zVertex[i]])
        for j in range(vertexDegree):
            eov1 = edgesOnVertex[i, j]
            eov2 = edgesOnVertex[i, (j+1)%vertexDegree]
            cov = cellsOnVertex[i, j]
            if cov >= 0 and eov1 >= 0 and eov2 >= 0:
                edge1 = np.asarray([xEdge[eov1], yEdge[eov1], zEdge[eov1]])
                edge2 = np.asarray([xEdge[eov2], yEdge[eov2], zEdge[eov2]])
                cell = np.asarray([xCell[cov], yCell[cov], zCell[cov]])
                kiteAreasOnVertex_new[i,j] = kiteAreasOnVertex_new[i,j] \
                                           + triangle_area(vtx, edge1, cell)
                kiteAreasOnVertex_new[i,j] = kiteAreasOnVertex_new[i,j] \
                                           + triangle_area(vtx, cell, edge2)
            else:
                kiteAreasOnVertex_new[i,j] = kiteAreasOnVertex[i,j] / scalefac**2

    mesh.variables['kiteAreasOnVertex'][:] = kiteAreasOnVertex_new
    kiteAreasOnVertex = kiteAreasOnVertex_new

    #
    # Create kiteAreasOnCell array
    #
    kiteAreasOnCell = np.zeros((nCells, maxEdges))

    for i in range(nCells):
        for j in range(nEdgesOnCell[i]):
            ii = verticesOnCell[i, j]
            k = find_neighbor(cellsOnVertex[ii, :], vertexDegree, i)
            if k >= 0:
                kiteAreasOnCell[i, j] = kiteAreasOnVertex[ii, k]
            else:
                print('Error finding cell %i in cellsOnVertex list for vertex %i (1-based)\n', i+1, ii+1)

    #
    # Normalize kiteAreasOnCell by areaCell
    #
    for i in range(nCells):
        v = 1.0 / areaCell[i]
        for j in range(nEdgesOnCell[i]):
            kiteAreasOnCell[i, j] = kiteAreasOnCell[i, j] * v

    #
    # weightsOnEdge
    #
    print('Recomputing weightsOnEdge...')
    weightsOnEdge_new = np.zeros((nEdges, maxEdges2))
    weightsOnEdge_new[:] = weightsOnEdge[:]
    for i in range(nEdges):
        j = 0

        if np.min(cellsOnEdge[i, :]) < 0:
            continue

        # Handle weights from first cell on this edge
        edge_id = i
        cell1 = cellsOnEdge[edge_id, 0]
        sum_r = 0.0
        rotated_kiteAreasOnCell = rotate_kiteAreasOnCell(kiteAreasOnCell[cell1, :],
                                                         edge_id,
                                                         edgesOnCell[cell1, :],
                                                         nEdgesOnCell[cell1],
                                                        )

        for ii in range(nEdgesOnCell[cell1] - 1):
            sum_r = sum_r + rotated_kiteAreasOnCell[ii]
            edge_id = edgesOnEdge[i, j]
            weightsOnEdge_new[i, j] = sign_for_edge(cell1, edge_id, cellsOnEdge) \
                                    * (0.5 - sum_r) * dvEdge[edge_id] / dcEdge[i]
            j = j + 1

        # Handle weights from second cell on this edge
        edge_id = i
        cell2 = cellsOnEdge[edge_id, 1]
        sum_r = 0.0
        rotated_kiteAreasOnCell = rotate_kiteAreasOnCell(kiteAreasOnCell[cell2, :],
                                                         edge_id,
                                                         edgesOnCell[cell2, :],
                                                         nEdgesOnCell[cell2],
                                                        )

        for ii in range(nEdgesOnCell[cell2] - 1):
            sum_r = sum_r + rotated_kiteAreasOnCell[ii];
            edge_id = edgesOnEdge[i, j]
            weightsOnEdge_new[i, j] = -1.0 * sign_for_edge(cell2, edge_id, cellsOnEdge) \
                                    * (0.5 - sum_r) * dvEdge[edge_id] / dcEdge[i]
            j = j + 1

    mesh.variables['weightsOnEdge'][:] = weightsOnEdge_new
    weightsOnEdge = weightsOnEdge_new

    mesh.close()
