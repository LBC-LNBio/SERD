import os
import pathlib
from typing import Union, Optional, Literal, List, Dict
import numpy
import networkx
from pyKVFinder import read_vdw, read_xyz
from pyKVFinder.grid import _get_sincos, _get_dimensions
from scipy.spatial.distance import cdist

__all__ = [
    "read_vdw",
    "read_xyz",
    "read_pdb",
    "get_vertices",
    "_get_sincos",
    "_get_dimensions",
    "surface",
    "interface",
    "detect",
    "save",
    "save_session",
    "r2g",
    "g2pdb",
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


def _process_pdb_line(
    line: str, vdw: Dict[str, Dict[str, float]]
) -> List[Union[str, float, int]]:
    """Extracts ATOM and HETATM information of PDB line.

    Parameters
    ----------
    line : str
        A line of a valid PDB file
    vdw : Dict[str, Dict[str, Dict[str, float]]]
        A dictionary containing radii values.

    Returns
    -------
    atomic : List[Union[str, float, int]]
        A list with resnum, chain, resname, atom name, xyz coordinates and radius.
    """
    # Get PDB infomation
    atom = line[12:16].strip()
    resname = line[17:20].strip()
    resnum = line[22:27].strip()
    chain = line[21]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    atom_symbol = line[76:78].strip().upper()

    # Get atom and radius from vdw
    if resname in vdw.keys() and atom in vdw[resname].keys():
        radius = vdw[resname][atom]
    else:
        radius = vdw["GEN"][atom_symbol]

    # Prepare output
    atomic = [resnum, chain, resname, atom, x, y, z, radius]

    return atomic


def read_pdb(
    fn: Union[str, pathlib.Path], vdw: Optional[Dict[str, Dict[str, float]]] = None
) -> numpy.ndarray:
    """Reads PDB file into numpy.ndarrays.

    Parameters
    ----------
    fn : Union[str, pathlib.Path]
        A path to PDB file.
    vdw : Dict[str, Dict[str, float]], optional
        A dictionary containing radii values, by default None. If None, use output of `pyKVFinder.read_vdw()`.

    Returns
    -------
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.

    Raises
    ------
    TypeError
        `fn` must be a string or a pathlib.Path.

    Note
    ----
    The van der Waals radii file defines the radius values for each atom
    by residue and when not defined, it uses a generic value based on the
    atom type. The function by default loads the built-in van der Waals radii
    file: `vdw.dat`.
    """
    # Check arguments
    if type(fn) not in [str, pathlib.Path]:
        raise TypeError("`fn` must be a string or a pathlib.Path.")

    # Define default vdw file
    if vdw is None:
        vdw = read_vdw(None)

    # Create lists
    atomic = []

    with open(fn, "r") as f:
        for line in f.readlines():
            if line[:4] == "ATOM" or line[:6] == "HETATM":
                atomic.append(_process_pdb_line(line, vdw))

    return numpy.asarray(atomic)


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
    from pyKVFinder import get_vertices as gv

    # Get vertices
    vertices = gv(atomic, 2 * probe, 2 * step)

    return vertices


def surface(
    atomic: numpy.ndarray,
    surface_representation: Literal["VDW", "SES", "SAS"] = "SES",
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
    surface_representation : Literal["VDW", "SES", "SAS"], optional
        Surface representation. Keywords options are VDW (van der Waals surface), SES (Solvent Excluded Surface)
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
        `surface_representation` must be a `VDW`, `SES` or `SAS`.
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
    if surface_representation not in ["VDW", "SES", "SAS"]:
        raise TypeError("`surface_representation` must be a `VDW`, `SES` or `SAS`.")
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
    if surface_representation == "VDW":
        if verbose:
            print("> Surface representation: van der Waals (vdW).")
        probe = 0.0
        surface_representation = True
    else:
        if probe == 0.0:
            raise ValueError(
                f"`probe` must be a positive real number, when {surface_representation} is set."
            )
        if surface_representation == "SES":
            if verbose:
                print("> Surface representation: Solvent Excluded Surface (SES).")
            surface_representation = True
        elif surface_representation == "SAS":
            if verbose:
                print("> Surface representation: Solvent Accessible Surface (SAS).")
            surface_representation = False

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
        surface_representation,
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
        surface,
        atominfo,
        xyzr,
        vertices[0],
        sincos,
        step,
        probe + step / 2,
        nthreads,
        verbose,
    )

    # Process residues
    residues = _process_residues(residues)

    return residues


def detect(
    target: Union[str, pathlib.Path],
    surface_representation: Literal["VDW", "SES", "SAS"] = "SES",
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
    surface_representation : Literal["VDW", "SES", "SAS"], optional
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
        `surface_representation` must be a `VDW`, `SES` or `SAS`.
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
    # Check arguments types
    if type(target) not in [str, pathlib.Path]:
        raise TypeError("`target` must be a string or a pathlib.Path.")
    if surface_representation not in ["VDW", "SES", "SAS"]:
        raise TypeError("`surface_representation` must be a `VDW`, `SES` or `SAS`.")
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

    # Read van der Waals radii dictionary
    vdw = read_vdw(vdw)

    # Read target biomolecule
    if target.endswith(".pdb"):
        atomic = read_pdb(target, vdw)
    elif target.endswith(".xyz"):
        atomic = read_xyz(target, vdw)
    else:
        raise ValueError("`target` must be .pdb or .xyz.")

    # Define solvent-exposed surface
    solvsurf = surface(atomic, surface_representation, step, probe, nthreads, verbose)

    # Define solvent-exposed residues
    residues = interface(
        solvsurf, atomic, ignore_backbone, step, probe, nthreads, verbose
    )

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
        A path to a PyMOL session file, by default "residues.pse".
    """
    import warnings
    from pymol import cmd

    if len(residues) < 1:
        warnings.warn("`residues` is an empty list.")
        return

    # Prepare target name
    target_name = os.path.basename(os.path.normpath(target)).replace(".pdb", "")

    # Load target biomolecule
    cmd.reinitialize()
    cmd.load(target, target_name)
    cmd.hide("everything", target_name)
    cmd.show("sticks", target_name)

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
    cmd.color("red", "residues")
    cmd.disable(target_name)
    cmd.enable(target_name)

    # Save PyMOL session
    cmd.save(fn)


def _get_atom_selection(
    atomic: numpy.ndarray, selection: Literal["CA", "CB"]
) -> numpy.ndarray:
    """Select atoms in atomic data.

    Parameters
    ----------
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    selection : {"CA", "CB"}, optional
        Atomic selection, by default "CB". Keywords options are:

            * 'CA': Select alfa-carbon;

            * 'CB': Select beta-carbon, except for glycine which selects the alfa-carbon.

    Returns
    -------
    atomic
        A numpy array with selected atomic data.
    """
    if selection == "CB":
        mask = numpy.logical_or(
            (atomic[:, 2:4] == [["GLY", "CA"]]).all(-1), atomic[:, 3] == "CB"
        )
    else:
        mask = atomic[:, 3] == "CA"

    return atomic[mask]


def _calculate_distance(atomic: numpy.ndarray) -> numpy.ndarray:
    """Calculate distance between atoms.

    Parameters
    ----------
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.

    Returns
    -------
    distance : numpy.ndarray
        Distance between atoms of atomic based on indexes.
    """
    atomic = atomic.astype(float)
    distance = cdist(atomic, atomic, metric="euclidean")
    return distance


def _all_atoms_to_residues(
    atomic: numpy.ndarray, distance: numpy.ndarray
) -> numpy.ndarray:
    """Convert distance between all atoms to minimal distance between residues.

    Parameters
    ----------
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    distance : numpy.ndarray
        Distance between all atoms.

    Returns
    -------
    residues : numpy.ndarray
        Distance between residues.
    """
    _, i, c = numpy.unique(
        atomic[:, 0:3], return_index=True, return_counts=True, axis=0
    )
    counts = c[numpy.argsort(i)]
    idx = numpy.sort(i)
    n_res = len(idx)

    residues = numpy.array(
        [
            numpy.min(
                distance[idx[i] : idx[i] + counts[i], idx[j] : idx[j] + counts[j]]
            )
            for i in range(n_res)
            for j in range(n_res)
        ]
    ).reshape((n_res, n_res))

    return residues


def r2g(
    residues: List[List[str]],
    atomic: numpy.ndarray,
    selection: Literal["CA", "CB"] = "CB",
    cutoff: Optional[float] = None,
    intraresidual: bool = False,
    weighted_edges: bool = False,
) -> networkx.classes.graph.Graph:
    """Create a graph from a list of solvent-exposed residues.

    Parameters
    ----------
    residues : List[List[str]]
        A list of solvent-exposed residues.
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    selection : {"CA", "CB", "all"}, optional
        Atomic selection, by default "CB". Keywords options are:

            * 'CA': Select alfa-carbon;

            * 'CB': Select beta-carbon, except for glycine which selects the alfa-carbon;

            * 'all': Select all atoms, distance between residues are the smallest distance between the atoms of these residues.

    cutoff : Optional[float], optional
        A limit of distance to define an edge between two solvent-exposed residues, by default None.
        If None, cutoff depends on selection argument. If "CA", cutoff is 10.0. If "CB", cutoff is 8.0.
    intraresidual : bool, optional
        Whether to consider intraresidual contacts to create adjacency matrix, by default False.
    weighted_edges : bool, optional
        Whether to include the distances as weight of the edges.

    Returns
    -------
    networkx.classes.graph.Graph
        A graph of solvent-exposed residues with edges defined by a distance smaller than the cutoff.

    Raises
    ------
    ValueError
        `selection` must be `CA`, `CB`, or `all`.

    Note
    ----
    Cutoff for beta-carbon is based on CAPRI round 28. For more details, refer to
    https://www.ebi.ac.uk/msd-srv/capri/round28/round28.html.
    """
    # Check arguments
    if cutoff is None:
        if selection == "CA":
            cutoff = 10.0
        elif selection == "CB":
            cutoff = 8.0
        elif selection == "all":
            cutoff = 5.0
        else:
            raise ValueError("`selection` must be `CA`, `CB`, or `all`.")

    # Get atom selection
    if selection != "all":
        atomic = _get_atom_selection(atomic, selection=selection)

    # Keep solvent exposed residues
    # https://stackoverflow.com/questions/51352527/check-for-identical-rows-in-different-numpy-arrays
    atomic = atomic[(atomic[:, 0:3][:, None] == residues).all(-1).any(-1)]

    # Calculate distance
    distance = _calculate_distance(atomic[:, 4:7])

    if selection == "all":
        distance = _all_atoms_to_residues(atomic, distance)

    # Calculate adjacency matrix
    if intraresidual:
        adjacency = (distance < cutoff).astype(int)
    else:
        adjacency = numpy.logical_and(distance > 0.0, distance < cutoff).astype(int)

    # Create networkx.Graph
    G = networkx.Graph()
    G.add_edges_from(numpy.argwhere(adjacency))
    if weighted_edges:
        weighted_edges = [
            (edges[0], edges[1], distance[edges[0]][edges[1]])
            for edges in G.edges(data=True)
        ]
        G.add_weighted_edges_from(weighted_edges)
    else:
        weighted_edges = [(edges[0], edges[1], 1.0) for edges in G.edges(data=True)]
        G.add_weighted_edges_from(weighted_edges)

    # Add interresidual distances as edge weights
    if weighted_edges:
        weighted_edges = [
            (edges[0], edges[1], distance[edges[0]][edges[1]])
            for edges in G.edges(data=True)
        ]
        G.add_weighted_edges_from(weighted_edges)

    return G


def g2pdb(
    graph: networkx.classes.graph.Graph,
    atomic: numpy.ndarray,
    residues: List[List[str]],
    fn: Union[str, pathlib.Path] = "graph.pdb",
):
    """Save a graph to a PDB-formatted file. Each node are represented by the
    CA atom of the residue and edges are represented by CONECT record.

    Parameters
    ----------
    graph : networkx.classes.graph.Graph
        A graph of solvent-exposed residues with edges defined by a distance smaller than the cutoff.
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
        and radius) for each atom.
    residues : List[List[str]]
        A list of solvent-exposed residues.
    fn : Union[str, pathlib.Path], optional
        A path to a PDB file, by default "graph.pdb".
    """
    # Get atom information of solvent-exposed residues
    atomic = _get_atom_selection(atomic, selection="CA")
    atomic = atomic[(atomic[:, 0:3][:, None] == residues).all(-1).any(-1)]

    # Write nodes to pdb
    with open(fn, "w") as f:
        for n_atom, atom in enumerate(atomic):
            if n_atom in list(graph.nodes):
                atomname = atom[3].center(4)
                resname = atom[2].ljust(3)
                chain = atom[1].rjust(1)
                resnum = atom[0].rjust(4)
                x = str("%8.3f" % (float(atom[4]))).rjust(8)
                y = str("%8.3f" % (float(atom[5]))).rjust(8)
                z = str("%8.3f" % (float(atom[6]))).rjust(8)
                f.write(
                    f"ATOM  {n_atom+1:5d} {atomname} {resname} {chain}{resnum}    {x}{y}{z}  1.00100.00\n"
                )
        for edge in graph.edges:
            if edge[0] != edge[1]:
                f.write(f"CONECT{str(edge[0]+1).rjust(5)}{str(edge[1]+1).rjust(5)}\n")
        f.write("END\n")
