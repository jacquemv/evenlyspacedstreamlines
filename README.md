
![Illustration of streamlines](https://github.com/jacquemv/evenlyspacedstreamlines/blob/main/illustration.png?raw=true)

### Objective

The objective is to generate a set of evenly-spaced streamlines tangent to a given vector field (or orientation field) on a smooth triangulated surface. The vector field is assumed to be constant over each triangle. A streamline is defined here as a polygonal line on the triangulated surface such that each segment of that polygonal line lies on a triangle and is parallel to the vector associated with that triangle. The algorithm distributes streamlines over the surface in such a way that the minimal distance between streamlines never becomes smaller than a given radius $r$. The approach is inspired by Jobard and Lefer [3].

The structure of the algorithm, its performance and its limitations are discussed in our papers [1, 2], with applications to the visualization of fiber orientation in the human atria.

### Minimal example

Here is an example of streamline generation on a triangular mesh composed of 3 vertices and 1 triangle, so there is 1 orientation vector.
```python
from evenlyspacedstreamlines import evenly_spaced_streamlines
vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
triangles = [[0, 1, 2]]
orientation = [[1.0, 1.0, 0.0]]
radius = 0.4
list_of_lines, list_of_indices, infos = \
    evenly_spaced_streamlines(vertices, triangles, orientation, radius)
```
The output gives the xyz coordinates of 3 streamlines, the indices of the triangles the streamlines pass through (here always triangle 0), and some geometrical information such as streamline lengths:
```python
>>> list_of_lines
[array([[0.4803, 0.5197, 0.],
        [0.,     0.0394, 0.]]),
 array([[0.7674, 0.2326, 0.],
        [0.5348, 0.,     0.]]),
 array([[0.1963, 0.8037, 0.],
        [0.,     0.6074, 0.]])]
>>> list_of_indices
[array([0], dtype=int32), array([0], dtype=int32), array([0], dtype=int32)]
>>> infos
StreamlinesInfos(lengths=array([0.6793, 0.3289, 0.2776]), min_altitude=0.7071,
max_base=1.4142, neighborhood_size=0.0, euler_characteristic=1, 
random_seed=1698460061)
```


### Syntax

```python
list_of_lines, list_of_indices, infos = \
    evenly_spaced_streamlines(vertices, triangles, orientation, radius,
        orthogonal=False, oriented_streamlines=False,
        seed_points=32, seed_region=None, max_length=0,
        avoid_u_turns=True, max_angle=90, singularity_mask_radius=0.1,
        allow_tweaking_orientation=True,
        random_seed=0, parallel=True, num_threads=-1)
```

### Positional arguments

The arguments **vertices** and **triangles** define the triangulated surface, and **orientation** the vector field on that surface. The parameter **radius** specifies streamlines spacing.
- **vertices** ($n_v$-by-3 float array): x, y, z coordinates of the $n_v$ vertices
- **triangles** ($n_t$-by-3 int array): indices of the vertices of the $n_t$ triangles
- **orientation** ($n_t$-by-3 float array): orientation vector in each triangle
- **radius** (float): the distance between streamlines will be larger than the radius and smaller than 2*radius

### Optional keyword arguments

- **seed_points** (int): number of seed points tested to generate each streamline; the longest streamline is kept (default: 32)
- **seed_region** (int array): list of triangle indices among which seed points are picked (default: None, which means that all triangles are considered)
- **orthogonal** (bool): if True, rotate the orientation by 90 degrees (default: False)
- **oriented_streamlines** (bool): if True, streamlines only follow the vector field in the direction it points to (and not the opposite direction); the outputted streamlines are then oriented according to the vector field. If False (default), the orientation field is defined modulo pi instead of 2pi
- **allow_tweaking_orientation** (bool): if an orientation vector is parallel to an edge of the triangle, a small random perturbation is applied to that vector to satisfy the requirement (default: True); otherwise an exception is raised when the requirement is not satisfied
- **singularity_mask_radius** (float): when the orientation field has a singularity (e.g. focus or node), prevent streamlines from entering a sphere of radius 'singularity_mask_radius' x 'radius' around that singularity (default: 0.1)
- **max_length** (int): maximal number of iterations when tracing streamlines; it is needed because of nearly-periodic streamlines (default: 0, which means equal to the number of triangles)
- **max_angle** (float): stop streamline integration if the angle between two consecutive segments is larger than max_angle in degrees; 0 means straight line and 180 means U-turn (default: 90)
- **avoid_u_turns** (bool): restrict high curvatures by maintaining a lateral (perpendicular) distance of at least 'radius' between a segment of the streamline and the other segments of the same streamline; this automatically sets 'max_angle' to at most 90 degrees (default: True)
- **random_seed** (int): initialize the seed for pseudo-random number generation (default: seed based on clock)
- **parallel** (bool): if True (default), use multithreading wherever implemented
- **num_threads** (int): if possible, use that number of threads for parallel computing (default: let OpenMP choose)

### Outputs

The function ``evenly_spaced_streamlines`` returns a 3-tuple:
- **list_of_lines** (list of $n$-by-3 matrices): xyz coordinates of each of the streamlines generated
- **list_of_indices** (list of ($n$-1)-arrays): vectors of indices indicating for each line segment of the streamline in which triangle they lie
- **infos** (namedtuple): information about streamline generation with the following attributes:
    - 'lengths' ($n$-array of float): length of each streamline;
    - 'min_altitude' (float): minimum altitude over all triangles;
    - 'max_base' (float): maximum length of the base edge over all triangles;
    - 'neighborhood_size' (float): average number of triangles at a distance < 'radius' from any triangle
    - 'euler_characteristic' (int): = #vertices - #edges + #triangles
    - 'random_seed' (int): random seed used for random number generation

### Visualization

A simple way to visualize the output is:
```python
import matplotlib.pyplot as plt
plt.figure().add_subplot(projection='3d')
for line in list_of_lines:
    plt.plot(line[:, 0], line[:, 1], line[:, 2])
plt.show()
```
For more sophisticated visualizations, the module also provides a function ``streamlines_to_tubes`` to represent the streamlines as a set of tubes (described by triangulated surfaces).

### Implementation

The code is implemented in C++, interfaced and compiled using cython, and with a wrapper for python (``evenlyspacedstreamlines/wrapper.py``). It was designed for meshes with a number of vertices of the order of 100k and a radius $r$ not too large as compared to mesh edge length.

### Installation

The package can be installed using the command ``pip install evenlyspacedstreamlines`` (on Windows, a compiler such as Microsoft Visual C++ is required).

If the code is downloaded from github, local installation on Linux is done by running ``make local`` and including the directory 'evenlyspacedstreamlines' in the PYTHONPATH environment variable. The easiest way to check if the installation worked is:
```python
from evenlyspacedstreamlines import test; test()
```
Tested using Anaconda 2023.09 (python 3.11) on Linux and Windows.

### Acknowledgements

This work was supported by the Natural Sciences and Engineering Research
Council of Canada (NSERC grant RGPIN-2020-05252).

### References

1. V. Jacquemet. Improved algorithm for generating evenly-spaced streamlines on a triangulated surface (*in preparation*), 2023.

2. A. Saliani, A. Tsikhanovich, V. Jacquemet. [Visualization of interpolated atrial fiber orientation using evenly-spaced streamlines](https://doi.org/10.1016/j.compbiomed.2019.103349), *Comput. Biol. Med.* 2019, vol. 11, pp. 103349. 

3. B. Jobard, W. Lefer. [Creating evenly-spaced streamlines of arbitrary density](https://link.springer.com/chapter/10.1007/978-3-7091-6876-9_5). In Visualization in Scientific Computing’97, pp. 43–55. Springer, 1997.