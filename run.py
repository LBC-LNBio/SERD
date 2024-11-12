# %%
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

    def calculate_atom_depth(self) -> pandas.DataFrame:
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
        data = numpy.c_[self.atomic[:, 0:4], atom_depth]

        return pandas.DataFrame(
            data, columns=["ResidueNumber", "Chain", "ResidueName", "AtomName", "AtomicDepth"], index=numpy.arange(1, len(data) + 1)
        )

    def get_interface(self):
        return interface(self.surface)



if __name__ == "__main__":

    structure = Structure()
    structure.load("examples/1FMO.pdb")
    structure.model_surface(type="SES", step=0.6, probe=1.4)
    atom_depth = structure.calculate_atom_depth()
    print(atom_depth)

# if __name__ == "__main__":
#     import argparse
#     parser = argparse.ArgumentParser(description="SERD")
#     parser.add_argument("pdb", type=str, help="Path to PDB file")
#     parser.add_argument("--vdw", type=str, default=None, help="Path to VDW file (optional)")
#     args = parser.parse_args()
#     print(args)

#     structure = Structure(vdw=args.vdw)
#     structure.load(args.pdb)
#     structure.model_surface()
#     print(structure.get_interface())

# %%
