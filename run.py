from SERD import read_vdw, read_pdb, get_vertices, surface, interface, _get_sincos
import numpy
import pandas
from typing import Dict, List, Optional


class Surface(object):

    def __init__(
        self, grid: numpy.ndarray, step: float, probe: float, vertices: numpy.ndarray
    ):
        self.grid = grid
        self.step = step
        self.probe = probe
        self.vertices = vertices
        self.coordinates = self._get_coordinates(grid, step, vertices)

    def _get_coordinates(
        self, grid: numpy.ndarray, step: float, vertices: numpy.ndarray
    ) -> numpy.ndarray:
        """
        Convert the grid representation of the surface to 3D Cartesian coordinates.

        Parameters
        ----------
        grid : numpy.ndarray
            The grid representation of the surface.
        step : float
            The step size used to model the surface.
        vertices : numpy.ndarray
            The vertices of the bounding box. P1: origin, P2: x-axis, P3: y-axis, P4: z-axis.

        Returns
        -------
        numpy.ndarray
            The 3D Cartesian coordinates of the surface.
        """
        indexes = numpy.argwhere(grid == 1)

        # P1, P2, P3, P4 (origin, x-axis, y-axis, z-axis)
        P1, _, _, _ = vertices

        # Calculate sin and cos for each axis
        sincos = _get_sincos(vertices)

        # Convert grid to 3D Cartesian coordinates
        xaux, yaux, zaux = (indexes * step).T

        x = (
            (xaux * sincos[3])
            + (yaux * sincos[0] * sincos[2])
            - (zaux * sincos[1] * sincos[2])
            + P1[0]
        )
        y = (yaux * sincos[1]) + (zaux * sincos[0]) + P1[1]
        z = (
            (xaux * sincos[2])
            - (yaux * sincos[0] * sincos[3])
            + (zaux * sincos[1] * sincos[3])
            + P1[2]
        )

        # Prepare 3D coordinates
        coordinates = numpy.array([x, y, z]).T

        return coordinates


class Structure(object):

    def __init__(self, vdw: Optional[str] = None, **kwargs):
        self.__dict__.update(kwargs)
        self.vdw = read_vdw(vdw)
        self.atomic = None
        self.surface = None

    def load(self, path: str):
        """
        Load the atomic data from a PDB file.

        Parameters
        ----------
        path : str
            The path to the PDB file.
        """
        self.atomic = read_pdb(path)

    def model_surface(self, type: str = "SES", step: float = 0.6, probe: float = 1.4):
        """
        Model the surface of the structure using the atomic data.
        The surface is modeled using the Solvent Excluded Surface (SES) or Solvent Accessible Surface (SAS) method. The SES method is used by default.

        Parameters
        ----------
        type : str, optional
            The type of surface to model, either 'SES' or 'SAS', by default 'SES'.
            SES: Solvent Excluded Surface. SAS: Solvent Accessible Surface.
        step : float, optional
            The step size used to model the surface, by default 0.6.
        probe : float, optional
            The radius of the probe used to model the surface, by default 1.4.

        Raises
        ------
        ValueError
            If no atomic data is loaded, raise an error.
        """
        if self.atomic is None:
            raise ValueError("No atomic data loaded. Please run .load() first.")

        # Calculate vertices of the bounding box
        vertices = get_vertices(self.atomic)

        # Model surface representation
        _surface = surface(
            self.atomic, surface_representation=type, step=step, probe=probe
        )
        self.surface = Surface(_surface, step, probe, vertices)

    def _get_interface(self, ignore_backbone: bool = True) -> Dict[str, List[int]]:
        return interface(
            self.surface.grid,
            self.atomic,
            ignore_backbone=ignore_backbone,
            step=self.surface.step,
            probe=self.surface.probe,
        )

    def atom_depth(self) -> pandas.DataFrame:
        """
        Calculate the depth of each atom in the structure. The atom radius is subtracted from the minimum distance to the surface.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the depth of each atom in the structure.
        """
        if surface is None:
            raise ValueError(
                "No surface data loaded. Please run .model_surface() first."
            )

        # Get coordinates from atomic
        atomic_coordinates = self.atomic[:, 4:7].astype(float)

        # Calculate distances between surface and atomic coordinates
        distances = numpy.sqrt(
            (
                (
                    self.surface.coordinates[:, numpy.newaxis, :]
                    - atomic_coordinates[numpy.newaxis, :, :]
                )
                ** 2
            ).sum(axis=2)
        )

        # Get minimum distance for each atom
        atom_depth = distances.min(axis=0) - self.atomic[:, 7].astype(float)

        # Prepare data
        data = pandas.DataFrame(
            self.atomic[:, 0:4],
            columns=["ResidueNumber", "Chain", "ResidueName", "AtomName"],
        )
        data["AtomicDepth"] = atom_depth

        return data

    def residue_depth(
        self,
        metric: str = "minimum",
        keep_only_interface: bool = False,
        ignore_backbone: bool = True,
    ) -> pandas.DataFrame:
        """
        Calculate the depth of each residue in the structure. The residue depth is calculated as the minimum or centroid of the atoms in the residue.

        Parameters
        ----------
        metric : str, optional
            The metric used to calculate the residue depth, either 'minimum', 'centroid', by default 'minimum'.
        keep_only_interface : bool, optional
            Whether to keep only residues at the interface, by default False.
        ignore_backbone : bool, optional
            Whether to ignore backbone atoms for defining the interface, by default True.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the depth of each residue in the structure.
        """
        # Calculate atom depth
        atom_depth = self.atom_depth()

        # Keep only residues at the interface
        if keep_only_interface:
            interface_residues = pandas.DataFrame(
                self._get_interface(ignore_backbone=ignore_backbone),
                columns=["ResidueNumber", "Chain", "ResidueName"],
            )

            # Keep only the interface
            atom_depth = pandas.merge(
                atom_depth,
                interface_residues,
                on=["ResidueNumber", "Chain", "ResidueName"],
                how="inner",
            )

        # Calculate residue depth
        if metric == "minimum":
            residue_depth = (
                atom_depth.groupby(["ResidueNumber", "Chain"], sort=False)
                .agg({"AtomicDepth": "min"})
                .reset_index()
            )
        elif metric == "centroid":
            residue_depth = (
                atom_depth.groupby(["ResidueNumber", "Chain"], sort=False)
                .agg({"AtomicDepth": "mean"})
                .reset_index()
            )
        else:
            raise ValueError("Invalid metric. Please use 'minimum' or 'centroid'.")

        # Rename column for consistency
        residue_depth.rename(columns={"AtomicDepth": "ResidueDepth"}, inplace=True)

        return residue_depth


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="SERD")
    parser.add_argument("pdb", type=str, help="Path to PDB file")
    parser.add_argument(
        "--vdw", type=str, default=None, help="Path to VDW file (optional)"
    )
    parser.add_argument(
        "--step", type=float, default=0.6, help="Step size for surface modeling"
    )
    parser.add_argument(
        "--probe", type=float, default=1.4, help="Probe radius for surface modeling"
    )
    parser.add_argument(
        "--type",
        type=str,
        default="SES",
        choices=["SES", "SAS"],
        help="Type of surface to model (SES or SAS)",
    )
    parser.add_argument(
        "--keep_only_interface",
        action="store_true",
        help="Keep only residues at the interface",
    )
    parser.add_argument(
        "--ignore_backbone",
        action="store_true",
        help="Ignore backbone atoms for defining the interface",
    )
    parser.add_argument(
        "--metric",
        type=str,
        default="minimum",
        choices=["minimum", "centroid"],
        help="Metric used to calculate the residue depth",
    )
    parser.add_argument(
        "--no-file", action="store_true", help="Do not write output to file"
    )
    args = parser.parse_args()

    # Load structure
    structure = Structure(vdw=args.vdw)
    structure.load(args.pdb)

    # Model surface
    structure.model_surface(type="SES", step=0.6, probe=1.4)

    # Calculate atom depth
    atom_depth = structure.atom_depth()
    if not args.no_file:
        atom_depth.to_csv("atom_depth.csv")
    print(atom_depth)

    # Calculate residue depth
    residue_depth = structure.residue_depth(
        metric=args.metric,
        keep_only_interface=args.keep_only_interface,
        ignore_backbone=args.ignore_backbone,
    )
    print(residue_depth)
    if not args.no_file:
        residue_depth.to_csv("residue_depth.csv")
