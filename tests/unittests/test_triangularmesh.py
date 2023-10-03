import numpy as np

def create_reference_data():
    from mtlheartmesh.trisurf import TriSurf
    basename = '../data/test_mesh3d'
    S = TriSurf.read_off(basename+'.off')
    with open(basename+'.nbound_ref', 'wt') as file:
        file.write(f'{len(S.boundary_vertices)}')
    with open(basename+'.radius_ref', 'wt') as file:
        file.write(f'{S.face_altitudes.min():.9f}')
    C = S.face_vertex_connectivity
    with open(basename+'.inctri_ref', 'wt') as file:
        for i in range(C.shape[0]):
            file.write(str(C[i].nonzero()[1])[1:-1]+'\n')
    S.adjacent_faces.tofile(basename+'.adjtri_ref', sep=' ')
    B = np.unique(np.concatenate(S.boundary_vertices))
    B.tofile(basename+'.bound_ref', sep=' ')
    with open(basename+'.neightri_ref', 'wt') as file:
        for i in range(S.nt):
            L = np.unique(np.concatenate([C[k].nonzero()[1]
                                          for k in S.tri[i]]))
            file.write(str(L)[1:-1]+'\n')

def test_geometric_calculations():
    basename = 'test_mesh3d'
    basename_ref = '../data/test_mesh3d'

    # check radius
    h = np.fromfile(basename+'.radius', sep=' ')
    h_ref = np.fromfile(basename_ref+'.radius_ref', sep=' ')
    assert(abs(h - h_ref) < 1e-5)
    
    # check v_incident_tri
    with open(basename+'.inctri') as file:
        with open(basename_ref+'.inctri_ref') as file_ref:
            for i, (line, line_ref) in enumerate(zip(file.readlines(),
                                                     file_ref.readlines())):
                items = sorted([int(k) for k in line.split(' ')[:-1]])
                items_ref = sorted([int(k) 
                                    for k in line_ref.strip().split(' ') 
                                    if k])
                assert np.all(np.array(items) == items_ref)

    # check t_adjacent_tri
    A = np.fromfile(basename+'.adjtri', sep=' ').reshape((-1, 3))
    A_ref = np.fromfile(basename_ref+'.adjtri_ref', sep=' ').reshape((-1, 3))
    assert np.all(A.sort(1) == A_ref.sort(1))

    #check is_boundary_vertex
    B = np.fromfile(basename+'.bound', sep=' ').astype(int)
    B_ref = np.fromfile(basename_ref+'.bound_ref', sep=' ').astype(int)
    assert np.all(np.array(B) == B_ref)

    # check t_neigh_tri
    with open(basename+'.neightri', 'rt') as file:
        with open(basename_ref+'.neightri_ref') as file_ref:
            for i, (line, line_ref) in enumerate(zip(file.readlines(),
                                                     file_ref.readlines())):
                items = sorted([i] + [int(k) for k in line.split(' ')[:-1]])
                items_ref = sorted([int(k) 
                                    for k in line_ref.strip().split(' ') 
                                    if k])
                assert np.all(np.array(items) == items_ref)


#create_reference_data()
test_geometric_calculations()