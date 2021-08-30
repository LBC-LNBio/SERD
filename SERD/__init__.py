"""Solvent-exposed residues detection (SERD).

SESD detects solvent-exposed residues of a target biomolecule.

Python package
--------------
>>> import SERD

See also
--------
* GitHub repository: https://github.com/jvsguerra/SERD
"""

__name__ = "SERD"
__version__ = "0.1.0"
__license__ = "GNU GPL-3.0 License"

import os
import pathlib
import numpy
from pyKVFinder import read_vdw, read_pdb, read_xyz
from pyKVFinder.grid import _get_sincos, _get_dimensions
from typing import Union, Optional, Literal, List

__all__ = [
    "read_vdw",
    "read_pdb",
    "read_xyz",
    "get_vertices",
    "_get_sincos",
    "_get_dimensions",
    "surface",
    "interface",
    "detect",
]


def _process_residues(residues: List[str]) -> List[List[str]]:
    return [res.split("_") for res in list(dict.fromkeys(residues))]


def get_vertices(
    atomic: numpy.ndarray,
    probe: Union[float, int] = 1.4,
    step: Union[float, int] = 0.6,
) -> numpy.ndarray:
    """Gets 3D grid vertices.

    Parameters
    ----------
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    probe : Union[float, int], optional
        Probe size (A), by default 4.0.
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6.

    Returns
    -------
    vertices : numpy.ndarray
        A numpy.ndarray with xyz vertices coordinates
        (origin, X-axis, Y-axis, Z-axis).
    """
    from pyKVFinder import get_vertices

    # Get vertices
    vertices = get_vertices(atomic, 2 * probe, step)

    return vertices


def surface(
    atomic: numpy.ndarray,
    surface: Literal["VDW", "SES", "SAS"] = "SES",
    step: Union[float, int] = 0.6,
    probe: Union[float, int] = 1.4,
    nthreads: Optional[int] = None,
    verbose: bool = False,
):
    from _SERD import _surface

    # Check arguments types
    if type(atomic) not in [numpy.ndarray, list]:
        raise TypeError("`atomic` must be a list or a numpy.ndarray.")
    elif len(numpy.asarray(atomic).shape) != 2:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    elif numpy.asarray(atomic).shape[1] != 8:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    if surface not in ["VDW", "SES", "SAS"]:
        raise TypeError("`surface` must be a `VDW`, `SES` or `SAS`.")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a non-negative real number.")
    if type(probe) not in [float, int]:
        raise TypeError("`probe` must be a non-negative real number.")
    elif probe < 0.0:
        raise ValueError("`probe` must be a non-negative real number.")
    if nthreads is None:
        nthreads = os.cpu_count() - 1
    else:
        if type(nthreads) not in [int]:
            raise TypeError("`nthreads` must be a positive integer.")
        elif nthreads <= 0:
            raise ValueError("`nthreads` must be a positive integer.")
    if type(verbose) not in [bool]:
        raise TypeError("`verbose` must be a boolean.")

    # Convert types
    step = float(step) if type(step) is int else step
    probe = float(probe) if type(probe) is int else probe

    # If surface representation is the van der Waals surface, the probe must be 0.0
    if surface == "VDW":
        if verbose:
            print("> Surface representation: van der Waals (vdW).")
        probe = 0.0
        surface = True
    else:
        if probe == 0.0:
            raise ValueError(
                f"`probe` must be a positive real number, when {surface} is set."
            )
        if surface == "SES":
            if verbose:
                print("> Surface representation: Solvent Excluded Surface (SES).")
            surface = True
        elif surface == "SAS":
            if verbose:
                print("> Surface representation: Solvent Accessible Surface (SAS).")
            surface = False

    # Get vertices
    vertices = get_vertices(atomic, probe, step)

    # Get sincos
    sincos = _get_sincos(vertices)

    # Get dimensions
    nx, ny, nz = _get_dimensions(vertices, step)

    # Get size
    size = nx * ny * nz

    # Extract xyzr from atomic
    xyzr = atomic[:, 4:].astype(numpy.float64)

    # Identify solvent-exposed surface
    surface = _surface(
        size,
        nx,
        ny,
        nz,
        xyzr,
        vertices[0],
        sincos,
        step,
        probe,
        surface,
        nthreads,
        verbose,
    ).reshape(nx, ny, nz)

    return surface


def interface(
    surface: numpy.ndarray,
    atomic: numpy.ndarray,
    ignore_backbone: bool = True,
    step: Union[float, int] = 0.6,
    probe: Union[float, int] = 1.4,
    nthreads: Optional[int] = None,
    verbose: bool = False,
):
    from _SERD import _interface

    # Check arguments types
    if type(surface) not in [numpy.ndarray]:
        raise TypeError("`surface` must be a numpy.ndarray.")
    elif len(surface.shape) != 3:
        raise ValueError("`cavities` has the incorrect shape. It must be (nx, ny, nz).")
    if type(atomic) not in [numpy.ndarray, list]:
        raise TypeError("`atomic` must be a list or a numpy.ndarray.")
    elif len(numpy.asarray(atomic).shape) != 2:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    elif numpy.asarray(atomic).shape[1] != 8:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    if type(ignore_backbone) not in [bool]:
        raise TypeError("`ignore_backbone` must be a boolean.")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a non-negative real number.")
    if type(probe) not in [float, int]:
        raise TypeError("`probe` must be a non-negative real number.")
    elif probe < 0.0:
        raise ValueError("`probe` must be a non-negative real number.")
    if nthreads is None:
        nthreads = os.cpu_count() - 1
    else:
        if type(nthreads) not in [int]:
            raise TypeError("`nthreads` must be a positive integer.")
        elif nthreads <= 0:
            raise ValueError("`nthreads` must be a positive integer.")
    if type(verbose) not in [bool]:
        raise TypeError("`verbose` must be a boolean.")

    # Get vertices
    vertices = get_vertices(atomic, probe, step)

    # Get sincos
    sincos = _get_sincos(vertices)

    # Get dimensions
    nx, ny, nz = _get_dimensions(vertices, step)

    # Extract xyzr from atomic
    xyzr = atomic[:, 4:].astype(numpy.float64)

    # Extract atominfo from atomic
    atominfo = numpy.asarray(
        ([[f"{atom[0]}_{atom[1]}_{atom[2]}", atom[3]] for atom in atomic[:, :4]])
    )

    # Remove backbone from atominfo
    if ignore_backbone:
        mask = numpy.where(
            (atominfo[:, 1] != "C")
            & (atominfo[:, 1] != "CA")
            & (atominfo[:, 1] != "N")
            & (atominfo[:, 1] != "O")
        )
        atominfo = atominfo[
            mask[0],
        ]
        xyzr = xyzr[
            mask[0],
        ]

    # Prepare atominfo
    atominfo = atominfo[:, 0].tolist()

    # Detect solvent-exposed residues
    residues = _interface(
        surface, atominfo, xyzr, vertices[0], sincos, step, probe, nthreads, verbose
    )

    # Process residues
    residues = _process_residues(residues)

    return residues


def detect(
    target: Union[str, pathlib.Path],
    surface: Literal["VDW", "SES", "SAS"] = "SES",
    step: Union[float, int] = 0.6,
    probe: Union[float, int] = 1.4,
    vdw: Optional[Union[str, pathlib.Path]] = None,
    ignore_backbone: bool = True,
    nthreads: Optional[int] = None,
    verbose: bool = False,
):
    from _SERD import _detect

    # Check arguments types
    if type(target) not in [str, pathlib.Path]:
        raise TypeError("`target` must be a string or a pathlib.Path.")
    if surface not in ["VDW", "SES", "SAS"]:
        raise TypeError("`surface` must be a `VDW`, `SES` or `SAS`.")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a non-negative real number.")
    if type(probe) not in [float, int]:
        raise TypeError("`probe` must be a non-negative real number.")
    if vdw is not None:
        if type(vdw) not in [str, pathlib.Path]:
            raise TypeError("`vdw` must be a string or a pathlib.Path.")
    if nthreads is None:
        nthreads = os.cpu_count() - 1
    else:
        if type(nthreads) not in [int]:
            raise TypeError("`nthreads` must be a positive integer.")
        elif nthreads <= 0:
            raise ValueError("`nthreads` must be a positive integer.")
    if type(verbose) not in [bool]:
        raise TypeError("`verbose` must be a boolean.")

    # Convert types
    step = float(step) if type(step) is int else step
    probe = float(probe) if type(probe) is int else probe

    # If surface representation is the van der Waals surface, the probe must be 0.0
    if surface == "VDW":
        if verbose:
            print("> Surface representation: van der Waals (vdW).")
        probe = 0.0
        surface = True
    else:
        if probe == 0.0:
            raise ValueError(
                f"`probe` must be a positive real number, when {surface} is set."
            )
        if surface == "SES":
            if verbose:
                print("> Surface representation: Solvent Excluded Surface (SES).")
            surface = True
        elif surface == "SAS":
            if verbose:
                print("> Surface representation: Solvent Accessible Surface (SAS).")
            surface = False

    # Read van der Waals radii dictionary
    vdw = read_vdw(vdw)

    # Read target biomolecule
    if target.endswith(".pdb"):
        atomic = read_pdb(target, vdw)
    elif target.endswith(".xyz"):
        atomic = read_xyz(target, vdw)
    else:
        raise ValueError("`target` must be .pdb or .xyz.")

    # Get vertices
    vertices = get_vertices(atomic, probe, step)

    # Get sincos
    sincos = _get_sincos(vertices)

    # Get dimensions
    nx, ny, nz = _get_dimensions(vertices, step)

    # Get size
    size = nx * ny * nz

    # Extract xyzr from atomic
    xyzr = atomic[:, 4:].astype(numpy.float64)

    # Extract atominfo from atomic
    atominfo = numpy.asarray(
        ([[f"{atom[0]}_{atom[1]}_{atom[2]}", atom[3]] for atom in atomic[:, :4]])
    )

    # Remove backbone from atominfo
    if ignore_backbone:
        mask = numpy.where(
            (atominfo[:, 1] != "C")
            & (atominfo[:, 1] != "CA")
            & (atominfo[:, 1] != "N")
            & (atominfo[:, 1] != "O")
        )
        atominfo = atominfo[
            mask[0],
        ]
        xyzr = xyzr[
            mask[0],
        ]

    # Prepare atominfo
    atominfo = atominfo[:, 0].tolist()

    # Detect solvent-exposed residues
    residues = _detect(
        size,
        nx,
        ny,
        nz,
        atominfo,
        xyzr,
        vertices[0],
        sincos,
        step,
        probe,
        surface,
        nthreads,
        verbose,
    )

    # Process residues
    residues = _process_residues(residues)

    return residues
