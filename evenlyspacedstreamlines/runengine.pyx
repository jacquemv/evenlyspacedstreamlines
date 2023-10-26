#cython: language_level=3
import numpy as np

#-----------------------------------------------------------------------------
cdef extern from "engine.h":
    cdef cppclass StreamlineEngine nogil:
        StreamlineEngine()
        void deallocate()
        void initialize(int nv_, int nt_, double* ver_, int* tri_, 
                        double* orient_, int orthogonal, 
                        int allow_tweaking_orientation)
        void setup(double radius_, int max_length, int max_nb_seeds_, 
                   int avoid_u_turns_, double max_angle, 
                   int oriented_streamlines_,
                   double singularity_mask_radius,
                   unsigned int random_seed, int parallel_, int num_threads)
        void define_seed_region(int seed_region_size, int* seed_region)
        void run()

        int get_nb_streamlines()
        int get_streamline_size(int no)
        void get_streamline(int no, double* points, int* idx)
        void get_lengths(double* L)

        int error_code
        double min_altitude
        double max_base
        int euler_characteristic
        double neighborhood_mean_size
        unsigned int saved_random_seed

#-----------------------------------------------------------------------------
def run_engine(double[:, ::1] vertices, 
               int[:, ::1] faces, 
               double[:, ::1] orient,
               double radius, 
               int[::1] seed_region,
               int orthogonal=False,
               int oriented_streamlines=False,
               int seed_points=32, 
               int max_length=0,
               int avoid_u_turns=True,
               double max_angle=90,
               double singularity_mask_radius=0.1,
               int allow_tweaking_orientation=True,
               unsigned int random_seed=0,
               int parallel=True,
               int num_threads=-1):
    cdef:
        int nv = vertices.shape[0]
        int nt = faces.shape[0]
        StreamlineEngine engine
        int ns, i, j, n
        double[:, ::1] points_memview
        int[::1] indices_memview
        double[::1] lengths_memview
    
    engine.initialize(nv, nt, &vertices[0, 0], &faces[0, 0], &orient[0, 0],
                      orthogonal, allow_tweaking_orientation)
    if engine.error_code == 3:
        raise ValueError('Non-manifold surface: an edge of the triangulated '
                         'surface has >2 adjacent triangles.')
    
    engine.setup(radius, max_length, seed_points, avoid_u_turns, max_angle,
                 oriented_streamlines, singularity_mask_radius, random_seed, 
                 parallel, num_threads)
    if seed_region.size:
        engine.define_seed_region(seed_region.size, &seed_region[0])
    if engine.error_code == 1:
        raise ValueError('At least one triangle of the mesh has '
                         'an area of zero.')
    elif engine.error_code == 2:
        raise ValueError('At least one orientation vector is parallel '
                         'to an edge of its corresponding triangle or has '
                         'zero norm. Try with '
                         'allow_tweaking_orientation=True')

    with nogil:
        engine.run()

    ns = engine.get_nb_streamlines()
    list_of_lines = []
    list_of_indices = []
    lengths = np.empty(ns, dtype=np.float64)

    for i in range(ns):
        n = engine.get_streamline_size(i)
        points = np.empty((n, 3), dtype=np.float64)
        indices = np.empty(n-1, dtype=np.int32)
        points_memview = points
        indices_memview = indices
        engine.get_streamline(i, &points_memview[0, 0], &indices_memview[0])
        list_of_lines.append(points)
        list_of_indices.append(indices)
    
    lengths_memview = lengths
    if ns > 0:
        engine.get_lengths(&lengths_memview[0])

    infos = dict(
                lengths=lengths, 
                min_altitude=float(engine.min_altitude),
                max_base=float(engine.max_base),
                neighborhood_size=float(engine.neighborhood_mean_size),
                euler_characteristic=int(engine.euler_characteristic),
                random_seed=int(engine.saved_random_seed)
            )
    return list_of_lines, list_of_indices, infos

