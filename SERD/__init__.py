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
    "save",
    "save_session",
]


def _process_residues(residues: List[str]) -> List[List[str]]:
    """Process raw list of residues from _interface or _detect to a list of
    residue information (residue number, chain identifier and residue name).

    Parameters
    ----------
    residues : List[str]
        A list of residues with duplications and items separated by '_'.

    Returns
    -------
    List[List[str]]
       A list of residue information (residue number, chain identifier and
       residue name).
    """
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
    vertices = get_vertices(atomic, 2 * probe, 2 * step)

    return vertices


def surface(
    atomic: numpy.ndarray,
    surface: Literal["VDW", "SES", "SAS"] = "SES",
    step: Union[float, int] = 0.6,
    probe: Union[float, int] = 1.4,
    nthreads: Optional[int] = None,
    verbose: bool = False,
) -> numpy.ndarray:
    """Defines the solvent-exposed surface of a target biomolecule in a 3D grid.

    Parameters
    ----------
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    surface : Literal["VDW", "SES", "SAS"], optional
        Surface representation. Keywords options are VDW (van der Waals), SES (Solvent Excluded Surface)
        or SAS (Solvent Accessible Surface), by default "SES".
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6.
    probe : Union[float, int], optional
        Probe size (A) to define SES and SAS representations, by default 1.4.
    nthreads : Optional[int], optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    surface : numpy.ndarray
        Surface points in the 3D grid (surface[nx, ny, nz]).
        Surface array has integer labels in each positions, that are:

            * -1: solvent points;

            * 0: biomolecule points;

            * 1: solvent-exposed surface points.

            Enclosed regions are considered biomolecule points.

    Raises
    ------
    TypeError
        `atomic` must be a numpy.ndarray.
    ValueError
        `atomic` has incorrect shape. It must be (n, 8).
    TypeError
        `surface` must be a `VDW`, `SES` or `SAS`.
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.
    TypeError
        `probe` must be a non-negative real number.
    ValueError
        `probe` must be a non-negative real number.
    TypeError
        `nthreads` must be a positive integer.
    ValueError
        `nthreads` must be a positive integer.
    TypeError
        `verbose` must be a boolean.
    ValueError
        `probe` must be a positive real number, when SES or SAS is set.
    """
    from _SERD import _surface

    # Check arguments types
    if type(atomic) not in [numpy.ndarray]:
        raise TypeError("`atomic` must be a numpy.ndarray.")
    elif len(atomic.shape) != 2:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    elif atomic.shape[1] != 8:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    if surface not in ["VDW", "SES", "SAS"]:
        raise TypeError("`surface` must be a `VDW`, `SES` or `SAS`.")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")
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
) -> List[List[str]]:
    """Identifies the solvent-exposed residues based on a target solvent-exposed surface
    and atomic information of a biomolecule (residue number, chain identifier, residue
    name, xyz coordinates, radius).

    Parameters
    ----------
    surface : numpy.ndarray
        Surface points in the 3D grid (surface[nx, ny, nz]).
        Surface array has integer labels in each positions, that are:

            * -1: solvent points;

            * 0: biomolecule points;

            * 1: solvent-exposed surface points.

            Enclosed regions are considered biomolecule points.
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    ignore_backbone : bool, optional
        Whether to ignore backbone atoms (C, CA, N, O) when defining interface residues, by default True.
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6.
    probe : Union[float, int], optional
        Probe size (A) to define SES and SAS representations, by default 1.4.
    nthreads : Optional[int], optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    residues : List[List[str]]
        A list of solvent-exposed residues.

    Raises
    ------
    TypeError
        `surface` must be a numpy.ndarray.
    ValueError
        `surface` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `atomic` must be a numpy.ndarray.
    ValueError
        `atomic` has incorrect shape. It must be (n, 8).
    TypeError
        `ignore_backbone` must be a boolean.
    TypeError
        `step` must be a positive real number.
    TypeError
        `probe` must be a non-negative real number.
    ValueError
        `probe` must be a non-negative real number.
    TypeError
        `nthreads` must be a positive integer.
    ValueError
        `nthreads` must be a positive integer.
    TypeError
        `verbose` must be a boolean.
    """
    from _SERD import _interface

    # Check arguments types
    if type(surface) not in [numpy.ndarray]:
        raise TypeError("`surface` must be a numpy.ndarray.")
    elif len(surface.shape) != 3:
        raise ValueError("`surface` has the incorrect shape. It must be (nx, ny, nz).")
    if type(atomic) not in [numpy.ndarray]:
        raise TypeError("`atomic` must be a numpy.ndarray.")
    elif len(atomic.shape) != 2:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    elif atomic.shape[1] != 8:
        raise ValueError("`atomic` has incorrect shape. It must be (n, 8).")
    if type(ignore_backbone) not in [bool]:
        raise TypeError("`ignore_backbone` must be a boolean.")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")
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
    """Detect solvent-exposed residues of a target biomolecule.

    Parameters
    ----------
    target : Union[str, pathlib.Path]
        A path to PDB or XYZ file of a target biomolecular structure.
    surface : Literal["VDW", "SES", "SAS"], optional
        Surface representation. Keywords options are VDW (van der Waals), SES (Solvent Excluded Surface)
        or SAS (Solvent Accessible Surface), by default "SES".
    step : Union[float, int], optional
        Grid spacing (A), by default 0.6.
    probe : Union[float, int], optional
        Probe size (A) to define SES and SAS representations, by default 1.4.
    vdw : Optional[Union[str, pathlib.Path]], optional
        A path to a van der Waals radii file, by default None. If None, apply the built-in van der
        Waals radii file: `vdw.dat`.
    ignore_backbone : bool, optional
        Whether to ignore backbone atoms (C, CA, N, O) when defining interface residues, by default True.
    nthreads : Optional[int], optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    verbose : bool, optional
        Print extra information to standard output, by default False.

    Returns
    -------
    residues : List[List[str]]
        A list of solvent-exposed residues.

    Raises
    ------
    TypeError
        `target` must be a string or a pathlib.Path.
    TypeError
        `surface` must be a `VDW`, `SES` or `SAS`.
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.
    TypeError
        `probe` must be a non-negative real number.
    ValueError
        `probe` must be a non-negative real number.
    TypeError
        `vdw` must be a string or a pathlib.Path.
    TypeError
        `nthreads` must be a positive integer.
    ValueError
        `nthreads` must be a positive integer.
    TypeError
        `verbose` must be a boolean.
    ValueError
        `probe` must be a positive real number, when SES or SAS is set.
    ValueError
        `target` must be .pdb or .xyz.

    Note
    ----
    The van der Waals radii file defines the radius values for each
    atom by residue and when not defined, it uses a generic value
    based on the atom type (see pyKVFinder package).
    """
    from _SERD import _detect

    # Check arguments types
    if type(target) not in [str, pathlib.Path]:
        raise TypeError("`target` must be a string or a pathlib.Path.")
    if surface not in ["VDW", "SES", "SAS"]:
        raise TypeError("`surface` must be a `VDW`, `SES` or `SAS`.")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")
    if type(probe) not in [float, int]:
        raise TypeError("`probe` must be a non-negative real number.")
    elif probe < 0.0:
        raise ValueError("`probe` must be a non-negative real number.")
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


def save(residues: List[List[str]], fn: Union[str, pathlib.Path] = "residues.pickle"):
    """Save list of solvent-exposed residues to binary pickle file.

    Parameters
    ----------
    residues : List[List[str]]
        A list of solvent-exposed residues.
    fn : Union[str, pathlib.Path], optional
        A path to pickle file, by default "residues.pickle"
    """
    import pickle

    with open(fn, "wb") as f:
        pickle.dump(residues, f)


def save_session(
    target: Union[str, pathlib.Path],
    residues: List[List[str]],
    fn: Union[str, pathlib.Path] = "residues.pse",
):
    """Save a PyMOL session with the solvent-exposed residues (shown as red sticks) and
    the target biomolecular structure (shown as cartoon).

    Parameters
    ----------
    target : Union[str, pathlib.Path]
        A path to PDB or XYZ file of a target biomolecular structure.
    residues : List[List[str]]
        A list of solvent-exposed residues.
    fn : Union[str, pathlib.Path], optional
        A path to a PyMOL session file, by default "residues.pse"
    """
    import warnings
    from pymol import cmd

    if len(residues) < 1:
        warnings.warn("`residues` is an empty list.")
        return

    # Prepare target name
    target_name = os.path.basename(os.path.normpath("examples/1FMO.pdb")).replace(
        ".pdb", ""
    )

    # Load target biomolecule
    cmd.reinitialize()
    cmd.load(target, target_name)

    # Select residues
    command = f"{target_name} and"
    for residue in residues:
        resnum, chain, _ = residue
        command = f"{command} (resid {resnum} and chain {chain}) or"
    command = f"{command[:-3]}"
    cmd.select("res", command)

    # Create residues object
    cmd.create("residues", "res")
    cmd.delete("res")
    cmd.hide("everything", "residues")
    cmd.show("sticks", "residues")
    cmd.color("red", "residues")
    cmd.disable(target_name)
    cmd.enable(target_name)
    cmd.set("auto_zoom", 1)

    # Save PyMOL session
    cmd.save(fn)
