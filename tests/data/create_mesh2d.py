#### Not for testing. To be run once.
import numpy as np
from scipy.spatial import Delaunay
from mtlheartmesh.tri2d import mesh_domain
from mtlheartmesh.trisurf import TriSurf, manip
import matplotlib.pyplot as plt

def create_contour(scale=2):
    points = np.array([(0, 0), (10, 1), (11, 6), (2, 8)])
    contour = []
    for i in range(4):
        x1, x2 = points[i], points[(i+1) % 4]
        t = np.linspace(0, 1, round(scale*np.linalg.norm(x2-x1)))[:-1][:, None]
        contour.append((1-t)*x1 + t*x2)
    return np.vstack(contour)

def create_coarse_mesh():
    contour = create_contour()
    ver, tri = mesh_domain(contour, max_area=0.5)
    S = TriSurf((ver, tri))
    S.summary()
    S.write_off('test_mesh2d.off')
    np.tile([2, 1, 0], (S.nt, 1)).tofile('test_mesh2d.orient', sep=' ')

    R = S.face_center - np.array([[5, 0, 0]])
    orient = np.column_stack((R[:, 0], R[:, 1], R[:, 2])).round(3)
    orient.tofile('test_mesh2d.orient2', sep=' ')

def create_fine_mesh():
    contour = create_contour(50)
    ver, tri = mesh_domain(contour, max_area=0.0004)
    S = TriSurf((ver, tri))
    S.summary()
    S.write_off('test_finemesh2d.off')

    R = S.face_center - np.array([[5, 0, 0]])
    orient = np.column_stack((R[:, 0], R[:, 1], R[:, 2])).round(3)
    orient.tofile('test_finemesh2d.orient', sep=' ')
    plt.triplot(S.ver[:, 0], S.ver[:, 1], S.tri)
    plt.show()

def create_annulus():
    n = 20
    r1, r2 = 1, 1.3

    theta = np.linspace(0, 2*np.pi, n+1)[:-1]
    x1, y1 = r1*np.cos(theta), r1*np.sin(theta)
    theta += theta[1]/2
    x2, y2 = r2*np.cos(theta), r2*np.sin(theta)

    ver = np.vstack((
        np.column_stack((x1, y1, y1*0)),
        np.column_stack((x2, y2, y2*0))
    ))
    I = np.arange(n)
    tri = np.vstack(( 
        np.column_stack((I, (I+1)%n, I+n)),
        np.column_stack(((I+1)%n, (I+1)%n + n, I+n))
    )).astype(np.int32)
    TriSurf((ver, tri)).write_off('test_annulus.off')
    
    G = (ver[tri[:, 0]] + ver[tri[:, 1]] + ver[tri[:, 2]]) / 3
    orient = np.column_stack((G[:, 1]+1.1e-9, -G[:, 0]+1e-9, 0*G[:, 0]))
    orient.tofile('test_annulus.orient', sep=' ')
    
    I = 2 * G[:, 0] + G[:, 1] < 0
    orient[I] = -orient[I]
    orient.tofile('test_annulus.orient2', sep=' ')

def create_moebius_old():
    u = np.linspace(0, 2*np.pi, 100)
    v = np.linspace(-1, 1, 20)
    u, v = np.meshgrid(u, v)
    u = u.flatten()
    v = v.flatten()

    tp = 1+0.5*v*np.cos(u/2)
    x = tp*np.cos(u)
    y = tp*np.sin(u)
    z = 0.5*v*np.sin(u/2)
    ver = np.column_stack((x, y, z))

    tri = Delaunay(np.vstack([u, v]).T).simplices
    I0 = np.where(u == 0)[0]
    I1 = np.where(u == 2*np.pi)[0]
    for i0, i1 in zip(I0, I1[::-1]):
        tri[tri == i1] = i0
    tri[::2] = tri[::2, [0, 2, 1]]
    S = TriSurf((ver, tri))
    S = manip.remove_isolated_vertices(S)
    S.write_off('test_moebius.off')

    x, y = S.ver[:, 0], S.ver[:, 1]
    V = np.column_stack((y, -x, 0*x))
    orient = V[S.tri].mean(1)
    S.project_vector(orient)
    S.normalize(orient)
    orient.tofile('test_moebius.orient', sep=' ')

def create_moebius():
    u = np.linspace(0, 2*np.pi, 100)
    v = np.linspace(-1, 1, 20)
    u, v = np.meshgrid(u, v)

    tp = 1+0.5*v*np.cos(u/2)
    x = tp*np.cos(u)
    y = tp*np.sin(u)
    z = 0.5*v*np.sin(u/2)
    S = manip.structured_mesh_to_trisurf(x, y, z)
    S.write_off('test_moebius.off')

    x, y = S.ver[:, 0], S.ver[:, 1]
    V = np.column_stack((y, -x, 0*x))
    orient = V[S.tri].mean(1)
    S.project_vector(orient)
    S.normalize(orient)
    orient.tofile('test_moebius.orient', sep=' ')

def create_klein():
    u = np.linspace(0, 2*np.pi, 200)
    v = np.linspace(0, 2*np.pi, 100) + np.pi/2
    u, v = np.meshgrid(u, v)

    a = 6*np.cos(u)*(1+np.sin(u))
    b = 16*np.sin(u)
    r = 4*(1-np.cos(u)/2)
    x = a+r*np.cos(u)*np.cos(v)
    y = b+r*np.sin(u)*np.cos(v)
    I = u > np.pi
    x[I] = a[I]+r[I]*np.cos(v[I]+np.pi)
    y[I] = b[I]
    z = r*np.sin(v)
    S, J = manip.structured_mesh_to_trisurf(x, y, z, return_indices=True)
    S.write_off('test_klein.off')

    vx = np.roll(x, 1, axis=0) - x
    vy = np.roll(y, 1, axis=0) - y
    vz = np.roll(z, 1, axis=0) - z
    V = np.column_stack((vx.ravel()[J], vy.ravel()[J], vz.ravel()[J]))

    V = V[S.tri]
    X = V[:, :, None, :] * V[:, :, :, None]
    Xm = X[:, 0] + X[:, 1] + X[:, 2]
    orient = np.linalg.eigh(Xm)[1][:, :, -1]

    #orient = V[S.tri].mean(1)
    S.project_vector(orient)
    orient = S.orthogonal_vector(orient)
    S.normalize(orient)
    orient.tofile('test_klein.orient', sep=' ')


def create_sphere():
    u = np.linspace(0, 2*np.pi, 100)
    v = np.linspace(0, np.pi, 50)
    u, v = np.meshgrid(u, v)

    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)

    S = manip.structured_mesh_to_trisurf(x, y, z, reverse=True)
    S.write_tri('dummy.tri')
    S.summary()

create_klein()
#create_sphere()

#create_moebius()
#create_annulus()
#create_fine_mesh()
