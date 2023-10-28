from collections import namedtuple
import numpy as np
from .runengine import run_engine

__all__ = ['evenly_spaced_streamlines', 'streamlines_to_tubes', 'test']

StreamlinesInfos = namedtuple('StreamlinesInfos', ['lengths', 'min_altitude', 
                              'max_base', 'neighborhood_size', 
                              'euler_characteristic', 'random_seed'])

#-----------------------------------------------------------------------------
def evenly_spaced_streamlines(vertices, triangles, orientation, radius, *,
                              orthogonal=False, 
                              oriented_streamlines=False,
                              seed_points=32, seed_region=None,
                              max_length=0,
                              avoid_u_turns=True, max_angle=90,
                              singularity_mask_radius=0.1,
                              allow_tweaking_orientation=True,
                              random_seed=0, parallel=True, num_threads=-1):
    """Generate a set of evenly spaced streamlines on a triangulated surface

    Args:
        vertices (nv-by-3 float array): x, y, z coordinates of the nv vertices
        triangles (nt-by-3 int array): indices of the vertices of the nt 
            triangles
        orientation (nt-by-3 float array): orientation vector in each triangle
        radius (float): the distance between streamlines will be larger than 
            the radius and smaller than 2*radius
    
    Optional keyword-only args:
        seed_points (int): number of seed points tested to generate each 
            streamline; the longest streamline is kept (default: 32)
        seed_region (int array): list of triangle indices among which seed 
            points are picked (default: None, which means that all triangles
            are considered)
        orthogonal (bool): if True, rotate the orientation by 90 degrees
            (default: False)
        oriented_streamlines (bool): if True, streamlines only follow the 
            vector field in the direction it points to (and not the opposite 
            direction); the outputted streamlines are then oriented according 
            to the vector field.
            If False (default), the orientation field is defined modulo pi 
            instead of 2pi
        allow_tweaking_orientation (bool): if an orientation vector is parallel
            to an edge of the triangle, a small random perturbation is applied 
            to that vector to satisfy the requirement (default: True); 
            otherwise an exception is raised when the requirement is not 
            satisfied
        singularity_mask_radius (float): when the orientation field has a 
            singularity (e.g. focus or node), prevent streamlines from entering
            a sphere of radius 'singularity_mask_radius' x 'radius' around 
            that singularity (default: 0.1)
        max_length (int): maximal number of iterations when tracing
            streamlines; it is needed because of nearly-periodic streamlines
            (default: 0, which means equal to the number of triangles)
        max_angle (float): stop streamline integration if the angle between
            two consecutive segments is larger than max_angle in degrees;
            0 means straight line and 180 means U-turn (default: 90)
        avoid_u_turns (bool): restrict high curvatures by maintaining a 
            lateral (perpendicular) distance of at least 'radius' between
            a segment of the streamline and the other segments of the same
            streamline; this automatically sets 'max_angle' to at most 
            90 degrees (default: True)
        random_seed (int): initialize the seed for pseudo-random number 
            generation (default: seed based on clock)
        parallel (bool): if True (default), use multithreading wherever 
            implemented
        num_threads (int): if possible, use that number of threads for parallel
            computing (default: let OpenMP choose)
    
    Returns:
        streamlines (list of n-by-3 matrices): xyz coordinates of each
            of the streamlines generated
        indices (list of (n-1)-arrays): vectors of indices indicating for each
            line segment of the streamline in which triangle they lie
        infos (namedtuple): information about streamline generation
            'lengths' (n-array): length of each streamline;
            'min_altitude' (float): minimum altitude over all triangles;
            'max_base' (float): maximum length of the base edge over all 
                triangles;
            'neighborhood_size' (float): average number of triangles at 
                a distance < 'radius' from any triangle
            'euler_characteristic' (int): = #vertices - #edges + #triangles
            'random_seed' (int): random seed used for random number generation
    """
    
    # ensure correct data types
    vertices = np.ascontiguousarray(vertices, dtype=np.float64)
    triangles = np.ascontiguousarray(triangles, dtype=np.int32)
    orientation = np.ascontiguousarray(orientation, dtype=np.float64)

    # for 2D meshes
    if vertices.shape[1] == 2:
        vertices = np.column_stack((vertices, np.zeros(vertices.shape[0])))
    if orientation.shape[1] == 2:
        orientation = np.column_stack((orientation, 
                                       np.zeros(orientation.shape[0])))

    # check number of columns
    if vertices.shape[1] != 3 or triangles.shape[1] != 3 \
                              or orientation.shape[1] != 3:
        raise ValueError("Arguments 'vertices', 'faces' and 'orient' "
                         "must all have 3 columns.")

    # check array sizes
    if orientation.shape[0] != triangles.shape[0]:
        raise ValueError("Arguments 'triangles' and 'orient' must have "
                         "the same number of rows.")
    
    # check bounds
    if triangles.min() < 0 or triangles.max() >= vertices.shape[0]:
        raise ValueError("Vertex indices out of bounds "
                         f"(must be between 0 and {vertices.shape[0]-1})")
    
    # seed region
    if seed_region is None:
        seed_region = np.empty(0, dtype=np.int32)
    else:
        seed_region = np.ascontiguousarray(seed_region, dtype=np.int32)

    # call C extension
    list_of_lines, list_of_indices, infos = run_engine(
        vertices, triangles, orientation, radius, 
        seed_region=seed_region,
        orthogonal=orthogonal,
        oriented_streamlines=oriented_streamlines,
        seed_points=seed_points, 
        max_length=max_length,
        avoid_u_turns=avoid_u_turns, 
        max_angle=max_angle,
        singularity_mask_radius=singularity_mask_radius,
        allow_tweaking_orientation=allow_tweaking_orientation,
        random_seed=random_seed, 
        parallel=parallel,
        num_threads=num_threads
    )
    return list_of_lines, list_of_indices, StreamlinesInfos(**infos)

#-------------------------------------------------------------------------------
def streamlines_to_tubes(surf_vertices, surf_triangles, lines, triangles, 
                         radius=0.01, nb_points=5, tapering=None):
    """Create a triangulated surface representing a set of tubes following
    the streamlines

    Args:
        surf_vertices (nv-by-3 array): vertices of the triangulated surface
        surf_triangles (nt-by-3 int array): triangulation of the surface (must 
            have consistent orientation)
        lines, triangles (lists): output of evenly_spaced_streamlines
        radius (float): radius of each tube
        nb_points (int): number of points around the section of a tube
            (default: 5)
        tapering (float): length scale of gradual reduction in tube radius at
            the extremities of the curves (default: None)
    
    Returns:
        vertices (float array), triangles (int array): triangulated surface 
        combining all the tubes
    """
    # vector normal to the surface
    P = surf_vertices[surf_triangles[:, 0], :]
    u = surf_vertices[surf_triangles[:, 1], :] - P
    v = surf_vertices[surf_triangles[:, 2], :] - P
    surf_normals = np.cross(u, v)
    surf_normals /= np.linalg.norm(surf_normals, axis=1)[:, None]

    # section of the tube
    ang = np.arange(nb_points) * 2*np.pi/nb_points
    x, y = np.cos(ang), np.sin(ang)

    all_ver, all_tri = [], []
    for line, face in zip(lines, triangles):
        
        # normal to the surface along the streamline
        normal = surf_normals[face]
        normal = np.row_stack((
            normal[0],
            normal[1:] + normal[:-1],
            normal[-1]
        ))
        normal /= np.maximum(np.linalg.norm(normal, axis=1)[:, None], 1e-12)

        # tangent to the streamline
        tangent = np.row_stack((
            line[1] - line[0],
            line[2:] - line[:-2],
            line[-1] - line[-2]
        ))
        tangent /= np.maximum(np.linalg.norm(tangent, axis=1)[:, None], 1e-12)

        # second vector orthogonal to the streamline
        orthogonal = np.cross(normal, tangent)

        # vertices of the tube
        tube_ver = x[:, None, None] * normal + y[:, None, None] * orthogonal
        if tapering:
            s = np.cumsum(np.sqrt(np.sum((line[1:]-line[:-1])**2, 1)))
            s = np.concatenate(([0], s))
            r = radius * np.expm1(-s/tapering) * np.expm1(-(s[-1]-s)/tapering)
            tube_ver *= r[None, :, None]
        else:
            tube_ver *= radius

        tube_ver += line
        tube_ver = tube_ver.reshape((-1, 3))
        all_ver.append(tube_ver)

        # triangulation of the tube
        tube_tri = []
        n = line.shape[0]
        I = np.arange(n-1, dtype=int)
        for k in range(nb_points):
            l = (k+1) % nb_points
            tube_tri.append(np.column_stack((I+k*n, I+k*n+1, I+l*n)))
            tube_tri.append(np.column_stack((I+k*n+1, I+l*n+1, I+l*n)))
        all_tri.append(np.vstack(tube_tri))
    
    # merge all tubes
    ver = np.vstack([v for v in all_ver])
    shift = np.cumsum([0] + [v.shape[0] for v in all_ver])
    tri = np.vstack([t+shift[i] for i, t in enumerate(all_tri)])
    return ver, tri

#-------------------------------------------------------------------------------
def test():
    """Run a simple example to check if the compiled code can be called.
    It should generate around 13 streamlines.
    """
    ver = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    tri = [[0, 1, 2]]
    orient = [[1, 1, 0]]
    radius = 0.1
    lines = evenly_spaced_streamlines(ver, tri, orient, radius)[0]
    print(f'{len(lines)} streamlines generated.')