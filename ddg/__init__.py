"""Discrete Differential Geometry benchmark on the unit sphere.

A compact, validated implementation of the DDG / DEC pipeline:

    mesh -> topology (incidence matrices) -> geometry (mass / Laplacian)
         -> heat diffusion -> geodesic distance (heat method)
         -> error & spectral analysis -> performance benchmarking

See the module docstrings and ``docs/theory.md`` for the mathematics, and
``scripts/run_benchmark.py`` for the end-to-end driver.
"""

from .mesh import TriMesh, icosphere, icosphere_for_target
from .topology import build_d0, build_d1, verify_d_squared
from .geometry import (
    face_areas,
    face_normals,
    interior_angles,
    vertex_areas,
    mass_matrix,
    cotangent_laplacian,
    cotangent_weights_per_edge,
)
from .heat import HeatSolver, point_source, total_heat
from .geodesic import (
    HeatMethodSolver,
    heat_geodesic,
    face_gradient,
    integrated_divergence,
)
from .spectral import eigenpairs, sphere_reference_spectrum, degree_for_index
from .analysis import exact_sphere_distance, error_metrics, convergence_rate

__all__ = [
    "TriMesh", "icosphere", "icosphere_for_target",
    "build_d0", "build_d1", "verify_d_squared",
    "face_areas", "face_normals", "interior_angles", "vertex_areas",
    "mass_matrix", "cotangent_laplacian", "cotangent_weights_per_edge",
    "HeatSolver", "point_source", "total_heat",
    "HeatMethodSolver", "heat_geodesic", "face_gradient",
    "integrated_divergence",
    "eigenpairs", "sphere_reference_spectrum", "degree_for_index",
    "exact_sphere_distance", "error_metrics", "convergence_rate",
]

__version__ = "0.1.0"
