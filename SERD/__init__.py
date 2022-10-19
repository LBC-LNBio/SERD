"""Solvent-exposed residues detection (SERD).

SESD detects solvent-exposed residues of a target biomolecule.

Python package
--------------
>>> import SERD

See also
--------
* GitHub repository: https://github.com/jvsguerra/SERD
* Documentation: https://github.com/jvsguerra/SERD/#api-reference
"""

__name__ = "SERD"
__version__ = "0.1.2"
__license__ = "GNU GPL-3.0 License"

try:
    from .SERD import *
except SyntaxError:
    pass
except ModuleNotFoundError:
    pass
