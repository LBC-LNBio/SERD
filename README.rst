####
SERD
####

.. image:: https://img.shields.io/pypi/v/SERD
    :target: https://pypi.org/project/SERD/

.. image:: https://img.shields.io/pypi/pyversions/SERD
    :target: https://pypi.org/project/SERD/

A Python package to detect solvent-exposed residues of a target biomolecule.

* `GitHub repository <https://github.com/jvsguerra/SERD>`_

************
Installation
************

To install the latest release on `PyPI <https://pypi.org/project/SERD>`_, 
run:

::

  pip install SERD

Or to install the latest developmental version, run:

::

  git clone https://github.com/jvsguerra/SERD.git
  pip install SERD

********
Tutorial
********

In this tutorial, we will use SERD on a catalytic subunit of a cAMP-dependent protein kinase (cADK) to detect its solvent-expoosed residues.

First of all, import SERD package on Python:

.. code:: python

  >>> import SERD

The full workflow for detecting solvent-exposed residues can be run with **SERD.detect** function:

.. code:: python
  
  >>> residues = SERD.detect('examples/1FMO.pdb')
  >>> print(residues)
  ['13', 'E', 'GLU'], ['14', 'E', 'SER'], ['15', 'E', 'VAL'], ['16', 'E', 'LYS'], ['17', 'E', 'GLU'], ['18', 'E', 'PHE'], ['19', 'E', 'LEU'], ['20', 'E', 'ALA'], ['21', 'E', 'LYS'], ['22', 'E', 'ALA'], ['23', 'E', 'LYS'], ['24', 'E', 'GLU'], ['25', 'E', 'ASP'], ['26', 'E', 'PHE'], ['27', 'E', 'LEU'], ['28', 'E', 'LYS'], ['29', 'E', 'LYS'], ['30', 'E', 'TRP'], ['31', 'E', 'GLU'], ['32', 'E', 'THR'], ['33', 'E', 'PRO'], ['34', 'E', 'SER'], ['35', 'E', 'GLN'], ['36', 'E', 'ASN'], ['37', 'E', 'THR'], ['38', 'E', 'ALA'], ['39', 'E', 'GLN'], ['40', 'E', 'LEU'], ['41', 'E', 'ASP'], ['42', 'E', 'GLN'], ['43', 'E', 'PHE'], ['44', 'E', 'ASP'], ['45', 'E', 'ARG'], ['46', 'E', 'ILE'], ['47', 'E', 'LYS'], ['48', 'E', 'THR'], ['49', 'E', 'LEU'], ['51', 'E', 'THR'], ['53', 'E', 'SER'], ['54', 'E', 'PHE'], ['56', 'E', 'ARG'], ['57', 'E', 'VAL'], ['58', 'E', 'MET'], ['59', 'E', 'LEU'], ['60', 'E', 'VAL'], ['61', 'E', 'LYS'], ['62', 'E', 'HIS'], ['63', 'E', 'LYS'], ['64', 'E', 'GLU'], ['65', 'E', 'SER'], ['67', 'E', 'ASN'], ['68', 'E', 'HIS'], ['69', 'E', 'TYR'], ['70', 'E', 'ALA'], ['71', 'E', 'MET'], ['72', 'E', 'LYS'], ['73', 'E', 'ILE'], ['74', 'E', 'LEU'], ['75', 'E', 'ASP'], ['76', 'E', 'LYS'], ['77', 'E', 'GLN'], ['78', 'E', 'LYS'], ['79', 'E', 'VAL'], ['80', 'E', 'VAL'], ['81', 'E', 'LYS'], ['82', 'E', 'LEU'], ['83', 'E', 'LYS'], ['84', 'E', 'GLN'], ['85', 'E', 'ILE'], ['86', 'E', 'GLU'], ['87', 'E', 'HIS'], ['88', 'E', 'THR'], ['89', 'E', 'LEU'], ['90', 'E', 'ASN'], ['91', 'E', 'GLU'], ['92', 'E', 'LYS'], ['93', 'E', 'ARG'], ['94', 'E', 'ILE'], ['95', 'E', 'LEU'], ['96', 'E', 'GLN'], ['97', 'E', 'ALA'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['100', 'E', 'PHE'], ['101', 'E', 'PRO'], ['102', 'E', 'PHE'], ['103', 'E', 'LEU'], ['104', 'E', 'VAL'], ['105', 'E', 'LYS'], ['106', 'E', 'LEU'], ['107', 'E', 'GLU'], ['108', 'E', 'PHE'], ['109', 'E', 'SER'], ['110', 'E', 'PHE'], ['111', 'E', 'LYS'], ['112', 'E', 'ASP'], ['113', 'E', 'ASN'], ['114', 'E', 'SER'], ['115', 'E', 'ASN'], ['116', 'E', 'LEU'], ['117', 'E', 'TYR'], ['118', 'E', 'MET'], ['119', 'E', 'VAL'], ['120', 'E', 'MET'], ['121', 'E', 'GLU'], ['122', 'E', 'TYR'], ['123', 'E', 'VAL'], ['124', 'E', 'ALA'], ['127', 'E', 'GLU'], ['128', 'E', 'MET'], ['129', 'E', 'PHE'], ['130', 'E', 'SER'], ['131', 'E', 'HIS'], ['132', 'E', 'LEU'], ['133', 'E', 'ARG'], ['134', 'E', 'ARG'], ['135', 'E', 'ILE'], ['137', 'E', 'ARG'], ['138', 'E', 'PHE'], ['139', 'E', 'SER'], ['140', 'E', 'GLU'], ['141', 'E', 'PRO'], ['142', 'E', 'HIS'], ['143', 'E', 'ALA'], ['144', 'E', 'ARG'], ['145', 'E', 'PHE'], ['146', 'E', 'TYR'], ['147', 'E', 'ALA'], ['148', 'E', 'ALA'], ['149', 'E', 'GLN'], ['150', 'E', 'ILE'], ['151', 'E', 'VAL'], ['152', 'E', 'LEU'], ['153', 'E', 'THR'], ['154', 'E', 'PHE'], ['155', 'E', 'GLU'], ['156', 'E', 'TYR'], ['157', 'E', 'LEU'], ['158', 'E', 'HIS'], ['159', 'E', 'SER'], ['160', 'E', 'LEU'], ['161', 'E', 'ASP'], ['162', 'E', 'LEU'], ['163', 'E', 'ILE'], ['164', 'E', 'TYR'], ['165', 'E', 'ARG'], ['166', 'E', 'ASP'], ['167', 'E', 'LEU'], ['168', 'E', 'LYS'], ['169', 'E', 'PRO'], ['170', 'E', 'GLU'], ['171', 'E', 'ASN'], ['172', 'E', 'LEU'], ['173', 'E', 'LEU'], ['174', 'E', 'ILE'], ['175', 'E', 'ASP'], ['176', 'E', 'GLN'], ['177', 'E', 'GLN'], ['179', 'E', 'TYR'], ['180', 'E', 'ILE'], ['181', 'E', 'GLN'], ['182', 'E', 'VAL'], ['183', 'E', 'THR'], ['184', 'E', 'ASP'], ['185', 'E', 'PHE'], ['187', 'E', 'PHE'], ['188', 'E', 'ALA'], ['189', 'E', 'LYS'], ['190', 'E', 'ARG'], ['191', 'E', 'VAL'], ['192', 'E', 'LYS'], ['194', 'E', 'ARG'], ['195', 'E', 'THR'], ['196', 'E', 'TRP'], ['197', 'E', 'TPO'], ['198', 'E', 'LEU'], ['199', 'E', 'CYS'], ['201', 'E', 'THR'], ['202', 'E', 'PRO'], ['203', 'E', 'GLU'], ['204', 'E', 'TYR'], ['205', 'E', 'LEU'], ['206', 'E', 'ALA'], ['207', 'E', 'PRO'], ['208', 'E', 'GLU'], ['209', 'E', 'ILE'], ['210', 'E', 'ILE'], ['211', 'E', 'LEU'], ['212', 'E', 'SER'], ['213', 'E', 'LYS'], ['215', 'E', 'TYR'], ['216', 'E', 'ASN'], ['217', 'E', 'LYS'], ['218', 'E', 'ALA'], ['219', 'E', 'VAL'], ['220', 'E', 'ASP'], ['221', 'E', 'TRP'], ['222', 'E', 'TRP'], ['223', 'E', 'ALA'], ['224', 'E', 'LEU'], ['226', 'E', 'VAL'], ['227', 'E', 'LEU'], ['228', 'E', 'ILE'], ['229', 'E', 'TYR'], ['230', 'E', 'GLU'], ['231', 'E', 'MET'], ['232', 'E', 'ALA'], ['233', 'E', 'ALA'], ['235', 'E', 'TYR'], ['236', 'E', 'PRO'], ['237', 'E', 'PRO'], ['238', 'E', 'PHE'], ['239', 'E', 'PHE'], ['240', 'E', 'ALA'], ['241', 'E', 'ASP'], ['242', 'E', 'GLN'], ['243', 'E', 'PRO'], ['244', 'E', 'ILE'], ['245', 'E', 'GLN'], ['246', 'E', 'ILE'], ['247', 'E', 'TYR'], ['248', 'E', 'GLU'], ['249', 'E', 'LYS'], ['250', 'E', 'ILE'], ['251', 'E', 'VAL'], ['252', 'E', 'SER'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['256', 'E', 'ARG'], ['257', 'E', 'PHE'], ['258', 'E', 'PRO'], ['259', 'E', 'SER'], ['260', 'E', 'HIS'], ['261', 'E', 'PHE'], ['262', 'E', 'SER'], ['263', 'E', 'SER'], ['264', 'E', 'ASP'], ['265', 'E', 'LEU'], ['266', 'E', 'LYS'], ['267', 'E', 'ASP'], ['268', 'E', 'LEU'], ['269', 'E', 'LEU'], ['270', 'E', 'ARG'], ['271', 'E', 'ASN'], ['272', 'E', 'LEU'], ['273', 'E', 'LEU'], ['274', 'E', 'GLN'], ['275', 'E', 'VAL'], ['276', 'E', 'ASP'], ['277', 'E', 'LEU'], ['278', 'E', 'THR'], ['279', 'E', 'LYS'], ['280', 'E', 'ARG'], ['281', 'E', 'PHE'], ['283', 'E', 'ASN'], ['284', 'E', 'LEU'], ['285', 'E', 'LYS'], ['286', 'E', 'ASN'], ['288', 'E', 'VAL'], ['289', 'E', 'ASN'], ['290', 'E', 'ASP'], ['291', 'E', 'ILE'], ['292', 'E', 'LYS'], ['293', 'E', 'ASN'], ['294', 'E', 'HIS'], ['295', 'E', 'LYS'], ['296', 'E', 'TRP'], ['297', 'E', 'PHE'], ['298', 'E', 'ALA'], ['299', 'E', 'THR'], ['300', 'E', 'THR'], ['301', 'E', 'ASP'], ['302', 'E', 'TRP'], ['303', 'E', 'ILE'], ['304', 'E', 'ALA'], ['305', 'E', 'ILE'], ['306', 'E', 'TYR'], ['307', 'E', 'GLN'], ['308', 'E', 'ARG'], ['309', 'E', 'LYS'], ['310', 'E', 'VAL'], ['311', 'E', 'GLU'], ['312', 'E', 'ALA'], ['313', 'E', 'PRO'], ['314', 'E', 'PHE'], ['315', 'E', 'ILE'], ['316', 'E', 'PRO'], ['317', 'E', 'LYS'], ['318', 'E', 'PHE'], ['319', 'E', 'LYS'], ['321', 'E', 'PRO'], ['323', 'E', 'ASP'], ['324', 'E', 'THR'], ['325', 'E', 'SER'], ['326', 'E', 'ASN'], ['327', 'E', 'PHE'], ['328', 'E', 'ASP'], ['329', 'E', 'ASP'], ['330', 'E', 'TYR'], ['331', 'E', 'GLU'], ['332', 'E', 'GLU'], ['333', 'E', 'GLU'], ['334', 'E', 'GLU'], ['335', 'E', 'ILE'], ['336', 'E', 'ARG'], ['337', 'E', 'VAL'], ['338', 'E', 'SEP'], ['339', 'E', 'ILE'], ['340', 'E', 'ASN'], ['341', 'E', 'GLU'], ['342', 'E', 'LYS'], ['343', 'E', 'CYS'], ['345', 'E', 'LYS'], ['346', 'E', 'GLU'], ['347', 'E', 'PHE'], ['348', 'E', 'THR'], ['349', 'E', 'GLU'], ['350', 'E', 'PHE']]

Still, users can save the list of solvent-exposed residues to a pickle binary file or save a PyMOL session to visualize these residues with the target biomolecular structure.

.. code:: python
  
  >>> SERD.save(residues, fn = 'residues.pickle')
  >>> SERD.save_session('examples/1FMO.pdb', residues, fn = 'residues.pse')

Furthermore, users can create a graph of solvent exposed residues, with edges based on a cutoff distance between a selection of atoms.

.. code:: python

  >>> atomic = SERD.read_pdb('examples/1FMO.pdb')
  >>> G = SERD.r2g(residues, atomic, selection="CB", cutoff=8.0)
  >>> G.number_of_nodes()
  263
  >>> G.number_of_edges()
  929

If users prefer, instead of running **SERD.detect** function, users can apply the detection of solvent-exposed residues in a step-by-step fashion. Below, we briefly describe this procedure.

**SERD.read_vdw** takes a vdW radii file (.dat) and returns a dictionary contaning radii values for each atom of each residue.

.. code:: python
  
  >>> vdw = SERD.read_vdw()
  >>> print(vdw)
  {'ALA': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB1': 1.487, '1HB': 1.487, 'HB2': 1.487, '2HB': 1.487, 'HB3': 1.487, '3HB': 1.487, 'C': 1.908, 'O': 1.6612}, 'ARG': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'HD2': 1.387, '1HD': 1.387, '2HD': 1.387, 'HD3': 1.387, 'HD1': 1.387, 'NE': 1.75, 'HE': 0.6, 'CZ': 1.908, 'NH1': 1.75, 'HH11': 0.6, '1HH1': 0.6, 'HH12': 0.6, '2HH1': 0.6, 'NH2': 1.75, 'HH21': 0.6, '2HH2': 0.6, 'HH22': 0.6, '1HH2': 0.6, 'C': 1.908, 'O': 1.6612}, 'ASH': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'OD2': 1.721, 'HD2': 0.0001, 'C': 1.908, 'O': 1.6612}, 'ASN': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'ND2': 1.824, 'HD21': 0.6, '1HD2': 0.6, 'HD22': 0.6, '2HD2': 0.6, 'C': 1.908, 'O': 1.6612}, 'ASP': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'OD2': 1.6612, 'C': 1.908, 'O': 1.6612}, 'CYM': {'N': 1.824, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB3': 1.387, 'HB2': 1.387, 'SG': 2.0, 'C': 1.908, 'O': 1.6612}, 'CYS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, '2HB': 1.387, '1HB': 1.387, 'HB3': 1.387, 'HB1': 1.387, 'SG': 2.0, 'HG': 0.6, 'C': 1.908, 'O': 1.6612}, 'CYX': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, 'SG': 2.0, 'C': 1.908, 'O': 1.6612}, 'GLH': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'HG2': 1.487, 'HG3': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'OE2': 1.721, 'HE2': 0.0001, 'C': 1.908, 'O': 1.6612}, 'GLN': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'NE2': 1.824, 'HE21': 0.6, '1HE2': 0.6, 'HE22': 0.6, '2HE2': 0.6, 'C': 1.908, 'O': 1.6612}, 'GLU': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'OE2': 1.6612, 'C': 1.908, 'O': 1.6612}, 'GLY': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA2': 1.387, 'HA1': 1.387, '1HA': 1.387, '2HA': 1.387, 'HA3': 1.387, 'C': 1.908, 'O': 1.6612}, 'HID': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'HIE': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'HIP': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'ILE': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.487, 'CG2': 1.908, 'HG21': 1.487, '1HG2': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG23': 1.487, '3HG2': 1.487, 'CG1': 1.908, 'HG12': 1.487, '2HG1': 1.487, 'HG13': 1.487, 'HG11': 1.487, '1HG1': 1.487, 'CD1': 1.908, 'HD11': 1.487, '1HD1': 1.487, 'HD12': 1.487, '2HD1': 1.487, 'HD13': 1.487, '3HD1': 1.487, 'C': 1.908, 'O': 1.6612}, 'LEU': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG': 1.487, 'CD1': 1.908, 'HD11': 1.487, '1HD1': 1.487, 'HD12': 1.487, '2HD1': 1.487, 'HD13': 1.487, '3HD1': 1.487, 'CD2': 1.908, 'HD21': 1.487, '1HD2': 1.487, 'HD22': 1.487, '2HD2': 1.487, 'HD23': 1.487, '3HD2': 1.487, 'C': 1.908, 'O': 1.6612}, 'LYN': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'HG2': 1.487, 'HG3': 1.487, 'CD': 1.908, 'HD2': 1.487, 'HD3': 1.487, 'CE': 1.908, 'HE2': 1.1, 'HE3': 1.1, 'NZ': 1.824, 'HZ2': 0.6, 'HZ3': 0.6, 'C': 1.908, 'O': 1.6612}, 'LYS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'HD2': 1.487, '1HD': 1.487, '2HD': 1.487, 'HD3': 1.487, 'HD1': 1.487, 'CE': 1.908, 'HE2': 1.1, '2HE': 1.1, 'HE3': 1.1, '1HE': 1.1, 'HE1': 1.1, 'NZ': 1.824, 'HZ1': 0.6, '1HZ': 0.6, 'HZ2': 0.6, '2HZ': 0.6, 'HZ3': 0.6, '3HZ': 0.6, 'C': 1.908, 'O': 1.6612}, 'MET': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.387, '2HG': 1.387, 'HG3': 1.387, 'HG1': 1.387, '1HG': 1.387, 'SD': 2.0, 'CE': 1.908, 'HE1': 1.387, '1HE': 1.387, 'HE2': 1.387, '2HE': 1.387, 'HE3': 1.387, '3HE': 1.387, 'C': 1.908, 'O': 1.6612}, 'PHE': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'HZ': 1.459, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'C': 1.908, 'O': 1.6612}, 'PRO': {'N': 1.824, 'CD': 1.908, 'HD2': 1.387, '1HD': 1.387, '2HD': 1.387, 'HD3': 1.387, 'HD1': 1.387, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CA': 1.908, 'HA': 1.387, 'C': 1.908, 'O': 1.6612}, 'SER': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, '2HB': 1.387, '1HB': 1.387, 'HB3': 1.387, 'HB1': 1.387, 'OG': 1.721, 'HG': 0.0001, 'C': 1.908, 'O': 1.6612}, 'THR': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, '1HG2': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG23': 1.487, '3HG2': 1.487, 'OG1': 1.721, 'HG1': 0.0001, 'C': 1.908, 'O': 1.6612}, 'TRP': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.85, 'CD1': 2.0, 'HD1': 1.409, 'NE1': 1.75, 'HE1': 0.6, 'CE2': 1.85, 'CZ2': 1.908, 'HZ2': 1.459, 'CH2': 1.908, 'HH2': 1.459, 'CZ3': 1.908, 'HZ3': 1.459, 'CE3': 1.908, 'HE3': 1.459, 'CD2': 1.85, 'C': 1.908, 'O': 1.6612}, 'TYR': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'OH': 1.721, 'HH': 0.0001, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'C': 1.908, 'O': 1.6612}, 'VAL': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.487, 'CG1': 1.908, 'CG2': 1.908, 'HG11': 1.487, '1HG2': 1.487, '1HG1': 1.487, 'HG21': 1.487, 'HG12': 1.487, '2HG1': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG13': 1.487, '3HG2': 1.487, '3HG1': 1.487, 'HG23': 1.487, 'C': 1.908, 'O': 1.6612}, 'HIS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'PTR': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'OH': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'SEP': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, '1HB': 1.387, '2HB': 1.387, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'TPO': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, 'HG22': 1.487, 'HG23': 1.487, '1HG2': 1.487, '2HG2': 1.487, '3HG2': 1.487, 'OG1': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'H2D': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'Y1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'T1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, 'HG22': 1.487, 'HG23': 1.487, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'S1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'GEN': {'AC': 2.0, 'AG': 1.72, 'AL': 2.0, 'AM': 2.0, 'AR': 1.88, 'AS': 1.85, 'AT': 2.0, 'AU': 1.66, 'B': 2.0, 'BA': 2.0, 'BE': 2.0, 'BH': 2.0, 'BI': 2.0, 'BK': 2.0, 'BR': 1.85, 'C': 1.66, 'CA': 2.0, 'CD': 1.58, 'CE': 2.0, 'CF': 2.0, 'CL': 1.75, 'CM': 2.0, 'CO': 2.0, 'CR': 2.0, 'CS': 2.0, 'CU': 1.4, 'DB': 2.0, 'DS': 2.0, 'DY': 2.0, 'ER': 2.0, 'ES': 2.0, 'EU': 2.0, 'F': 1.47, 'FE': 2.0, 'FM': 2.0, 'FR': 2.0, 'GA': 1.87, 'GD': 2.0, 'GE': 2.0, 'H': 0.91, 'HE': 1.4, 'HF': 2.0, 'HG': 1.55, 'HO': 2.0, 'HS': 2.0, 'I': 1.98, 'IN': 1.93, 'IR': 2.0, 'K': 2.75, 'KR': 2.02, 'LA': 2.0, 'LI': 1.82, 'LR': 2.0, 'LU': 2.0, 'MD': 2.0, 'MG': 1.73, 'MN': 2.0, 'MO': 2.0, 'MT': 2.0, 'N': 1.97, 'NA': 2.27, 'NB': 2.0, 'ND': 2.0, 'NE': 1.54, 'NI': 1.63, 'NO': 2.0, 'NP': 2.0, 'O': 1.69, 'OS': 2.0, 'P': 2.1, 'PA': 2.0, 'PB': 2.02, 'PD': 1.63, 'PM': 2.0, 'PO': 2.0, 'PR': 2.0, 'PT': 1.72, 'PU': 2.0, 'RA': 2.0, 'RB': 2.0, 'RE': 2.0, 'RF': 2.0, 'RH': 2.0, 'RN': 2.0, 'RU': 2.0, 'S': 2.09, 'SB': 2.0, 'SC': 2.0, 'SE': 1.9, 'SG': 2.0, 'SI': 2.1, 'SM': 2.0, 'SN': 2.17, 'SR': 2.0, 'TA': 2.0, 'TB': 2.0, 'TC': 2.0, 'TE': 2.06, 'TH': 2.0, 'TI': 2.0, 'TL': 1.96, 'TM': 2.0, 'U': 1.86, 'V': 2.0, 'W': 2.0, 'XE': 2.16, 'Y': 2.0, 'YB': 2.0, 'ZN': 1.39, 'ZR': 2.0}}

.. note::
 
  This step is only necessary if you are reading a custom van der Waals radii file to use in **SERD.read_pdb**.

.. seealso::

  `pyKVFinder.read_vdw <https://lbc-lnbio.github.io/pyKVFinder/_api_reference/read_vdw.html>`_

**SERD.read_pdb** or **SERD.read_xyz** take a target .pdb or .xyz file and returns the atomic information (residue number, chain identifier, residue name, xyz coordinates, radius).

.. code:: python

  >>> atomic = SERD.read_pdb('examples/1FMO.pdb')
  >>> print(atomic)
  array([['13', 'E', 'GLU', ..., '-15.642', '-14.858', '1.824'],
     ['13', 'E', 'GLU', ..., '-14.62', '-15.897', '1.908'],
     ['13', 'E', 'GLU', ..., '-13.357', '-15.508', '1.908'],
     ...,
     ['350', 'E', 'PHE', ..., '18.878', '-9.885', '1.908'],
     ['350', 'E', 'PHE', ..., '17.624', '-9.558', '1.908'],
     ['350', 'E', 'PHE', ..., '19.234', '-13.442', '1.69']],
    dtype='<U32')

.. note::
  
  The function takes the built-in dictionary, when the vdw argument is not specified. If you wish to use a custom van der Waals radii file, you must read it with pyKVFinder.read_vdw as shown earlier and pass it as **SERD.read_pdb(pdb, vdw=vdw)**.

**SERD.get_vertices** takes atomic information of a biomolecule (residue number, chain identifier, residue name, xyz coordinates, radius), and probe (`probe`) and grid spacing (`step`) that will be applied in the detection, and returns a NumPy array with vertice coordinates (origin, X-axis, Y-axis, Z-axis) of the 3D grid.

.. code:: python

  >>> vertices = SERD.get_vertices(atomic)
  >>> print(vertices)
  [[-19.311 -31.525 -30.206]
   [ 39.588 -31.525 -30.206]
   [-19.311  42.846 -30.206]
   [-19.311 -31.525  26.752]]

**SERD.surface** takes atomic information of a biomolecule (residue number, chain identifier, residue name, xyz coordinates, radius), and defines the solvent-exposed surface of a target biomolecule in a 3D grid.

.. code:: python

  >>> surface = SERD.surface(atomic)
  >>> print(surface)
  array([[[-1, -1, -1, ..., -1, -1, -1],
      [-1, -1, -1, ..., -1, -1, -1],
      [-1, -1, -1, ..., -1, -1, -1],
      ...,
      [-1, -1, -1, ..., -1, -1, -1],
      [-1, -1, -1, ..., -1, -1, -1],
      [-1, -1, -1, ..., -1, -1, -1]],

     ...,

     [[-1, -1, -1, ..., -1, -1, -1],
      [-1, -1, -1, ..., -1, -1, -1],
      [-1, -1, -1, ..., -1, -1, -1],
      ...,
      [-1, -1, -1, ..., -1, -1, -1],
      [-1, -1, -1, ..., -1, -1, -1],
      [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)

**SERD.interface** takes a target solvent-exposed surface (3D grid) and atomic information of a biomolecule (residue number, chain identifier, residue name, xyz coordinates, radius), and identifies the solvent-exposed residues.

.. code:: python

  >>> residues = SERD.interface(surface, atomic)
  >>> print(residues)
  ['13', 'E', 'GLU'], ['14', 'E', 'SER'], ['15', 'E', 'VAL'], ['16', 'E', 'LYS'], ['17', 'E', 'GLU'], ['18', 'E', 'PHE'], ['19', 'E', 'LEU'], ['20', 'E', 'ALA'], ['21', 'E', 'LYS'], ['22', 'E', 'ALA'], ['23', 'E', 'LYS'], ['24', 'E', 'GLU'], ['25', 'E', 'ASP'], ['26', 'E', 'PHE'], ['27', 'E', 'LEU'], ['28', 'E', 'LYS'], ['29', 'E', 'LYS'], ['30', 'E', 'TRP'], ['31', 'E', 'GLU'], ['32', 'E', 'THR'], ['33', 'E', 'PRO'], ['34', 'E', 'SER'], ['35', 'E', 'GLN'], ['36', 'E', 'ASN'], ['37', 'E', 'THR'], ['38', 'E', 'ALA'], ['39', 'E', 'GLN'], ['40', 'E', 'LEU'], ['41', 'E', 'ASP'], ['42', 'E', 'GLN'], ['43', 'E', 'PHE'], ['44', 'E', 'ASP'], ['45', 'E', 'ARG'], ['46', 'E', 'ILE'], ['47', 'E', 'LYS'], ['48', 'E', 'THR'], ['49', 'E', 'LEU'], ['51', 'E', 'THR'], ['53', 'E', 'SER'], ['54', 'E', 'PHE'], ['56', 'E', 'ARG'], ['57', 'E', 'VAL'], ['58', 'E', 'MET'], ['59', 'E', 'LEU'], ['60', 'E', 'VAL'], ['61', 'E', 'LYS'], ['62', 'E', 'HIS'], ['63', 'E', 'LYS'], ['64', 'E', 'GLU'], ['65', 'E', 'SER'], ['67', 'E', 'ASN'], ['68', 'E', 'HIS'], ['69', 'E', 'TYR'], ['70', 'E', 'ALA'], ['71', 'E', 'MET'], ['72', 'E', 'LYS'], ['73', 'E', 'ILE'], ['74', 'E', 'LEU'], ['75', 'E', 'ASP'], ['76', 'E', 'LYS'], ['77', 'E', 'GLN'], ['78', 'E', 'LYS'], ['79', 'E', 'VAL'], ['80', 'E', 'VAL'], ['81', 'E', 'LYS'], ['82', 'E', 'LEU'], ['83', 'E', 'LYS'], ['84', 'E', 'GLN'], ['85', 'E', 'ILE'], ['86', 'E', 'GLU'], ['87', 'E', 'HIS'], ['88', 'E', 'THR'], ['89', 'E', 'LEU'], ['90', 'E', 'ASN'], ['91', 'E', 'GLU'], ['92', 'E', 'LYS'], ['93', 'E', 'ARG'], ['94', 'E', 'ILE'], ['95', 'E', 'LEU'], ['96', 'E', 'GLN'], ['97', 'E', 'ALA'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['100', 'E', 'PHE'], ['101', 'E', 'PRO'], ['102', 'E', 'PHE'], ['103', 'E', 'LEU'], ['104', 'E', 'VAL'], ['105', 'E', 'LYS'], ['106', 'E', 'LEU'], ['107', 'E', 'GLU'], ['108', 'E', 'PHE'], ['109', 'E', 'SER'], ['110', 'E', 'PHE'], ['111', 'E', 'LYS'], ['112', 'E', 'ASP'], ['113', 'E', 'ASN'], ['114', 'E', 'SER'], ['115', 'E', 'ASN'], ['116', 'E', 'LEU'], ['117', 'E', 'TYR'], ['118', 'E', 'MET'], ['119', 'E', 'VAL'], ['120', 'E', 'MET'], ['121', 'E', 'GLU'], ['122', 'E', 'TYR'], ['123', 'E', 'VAL'], ['124', 'E', 'ALA'], ['127', 'E', 'GLU'], ['128', 'E', 'MET'], ['129', 'E', 'PHE'], ['130', 'E', 'SER'], ['131', 'E', 'HIS'], ['132', 'E', 'LEU'], ['133', 'E', 'ARG'], ['134', 'E', 'ARG'], ['135', 'E', 'ILE'], ['137', 'E', 'ARG'], ['138', 'E', 'PHE'], ['139', 'E', 'SER'], ['140', 'E', 'GLU'], ['141', 'E', 'PRO'], ['142', 'E', 'HIS'], ['143', 'E', 'ALA'], ['144', 'E', 'ARG'], ['145', 'E', 'PHE'], ['146', 'E', 'TYR'], ['147', 'E', 'ALA'], ['148', 'E', 'ALA'], ['149', 'E', 'GLN'], ['150', 'E', 'ILE'], ['151', 'E', 'VAL'], ['152', 'E', 'LEU'], ['153', 'E', 'THR'], ['154', 'E', 'PHE'], ['155', 'E', 'GLU'], ['156', 'E', 'TYR'], ['157', 'E', 'LEU'], ['158', 'E', 'HIS'], ['159', 'E', 'SER'], ['160', 'E', 'LEU'], ['161', 'E', 'ASP'], ['162', 'E', 'LEU'], ['163', 'E', 'ILE'], ['164', 'E', 'TYR'], ['165', 'E', 'ARG'], ['166', 'E', 'ASP'], ['167', 'E', 'LEU'], ['168', 'E', 'LYS'], ['169', 'E', 'PRO'], ['170', 'E', 'GLU'], ['171', 'E', 'ASN'], ['172', 'E', 'LEU'], ['173', 'E', 'LEU'], ['174', 'E', 'ILE'], ['175', 'E', 'ASP'], ['176', 'E', 'GLN'], ['177', 'E', 'GLN'], ['179', 'E', 'TYR'], ['180', 'E', 'ILE'], ['181', 'E', 'GLN'], ['182', 'E', 'VAL'], ['183', 'E', 'THR'], ['184', 'E', 'ASP'], ['185', 'E', 'PHE'], ['187', 'E', 'PHE'], ['188', 'E', 'ALA'], ['189', 'E', 'LYS'], ['190', 'E', 'ARG'], ['191', 'E', 'VAL'], ['192', 'E', 'LYS'], ['194', 'E', 'ARG'], ['195', 'E', 'THR'], ['196', 'E', 'TRP'], ['197', 'E', 'TPO'], ['198', 'E', 'LEU'], ['199', 'E', 'CYS'], ['201', 'E', 'THR'], ['202', 'E', 'PRO'], ['203', 'E', 'GLU'], ['204', 'E', 'TYR'], ['205', 'E', 'LEU'], ['206', 'E', 'ALA'], ['207', 'E', 'PRO'], ['208', 'E', 'GLU'], ['209', 'E', 'ILE'], ['210', 'E', 'ILE'], ['211', 'E', 'LEU'], ['212', 'E', 'SER'], ['213', 'E', 'LYS'], ['215', 'E', 'TYR'], ['216', 'E', 'ASN'], ['217', 'E', 'LYS'], ['218', 'E', 'ALA'], ['219', 'E', 'VAL'], ['220', 'E', 'ASP'], ['221', 'E', 'TRP'], ['222', 'E', 'TRP'], ['223', 'E', 'ALA'], ['224', 'E', 'LEU'], ['226', 'E', 'VAL'], ['227', 'E', 'LEU'], ['228', 'E', 'ILE'], ['229', 'E', 'TYR'], ['230', 'E', 'GLU'], ['231', 'E', 'MET'], ['232', 'E', 'ALA'], ['233', 'E', 'ALA'], ['235', 'E', 'TYR'], ['236', 'E', 'PRO'], ['237', 'E', 'PRO'], ['238', 'E', 'PHE'], ['239', 'E', 'PHE'], ['240', 'E', 'ALA'], ['241', 'E', 'ASP'], ['242', 'E', 'GLN'], ['243', 'E', 'PRO'], ['244', 'E', 'ILE'], ['245', 'E', 'GLN'], ['246', 'E', 'ILE'], ['247', 'E', 'TYR'], ['248', 'E', 'GLU'], ['249', 'E', 'LYS'], ['250', 'E', 'ILE'], ['251', 'E', 'VAL'], ['252', 'E', 'SER'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['256', 'E', 'ARG'], ['257', 'E', 'PHE'], ['258', 'E', 'PRO'], ['259', 'E', 'SER'], ['260', 'E', 'HIS'], ['261', 'E', 'PHE'], ['262', 'E', 'SER'], ['263', 'E', 'SER'], ['264', 'E', 'ASP'], ['265', 'E', 'LEU'], ['266', 'E', 'LYS'], ['267', 'E', 'ASP'], ['268', 'E', 'LEU'], ['269', 'E', 'LEU'], ['270', 'E', 'ARG'], ['271', 'E', 'ASN'], ['272', 'E', 'LEU'], ['273', 'E', 'LEU'], ['274', 'E', 'GLN'], ['275', 'E', 'VAL'], ['276', 'E', 'ASP'], ['277', 'E', 'LEU'], ['278', 'E', 'THR'], ['279', 'E', 'LYS'], ['280', 'E', 'ARG'], ['281', 'E', 'PHE'], ['283', 'E', 'ASN'], ['284', 'E', 'LEU'], ['285', 'E', 'LYS'], ['286', 'E', 'ASN'], ['288', 'E', 'VAL'], ['289', 'E', 'ASN'], ['290', 'E', 'ASP'], ['291', 'E', 'ILE'], ['292', 'E', 'LYS'], ['293', 'E', 'ASN'], ['294', 'E', 'HIS'], ['295', 'E', 'LYS'], ['296', 'E', 'TRP'], ['297', 'E', 'PHE'], ['298', 'E', 'ALA'], ['299', 'E', 'THR'], ['300', 'E', 'THR'], ['301', 'E', 'ASP'], ['302', 'E', 'TRP'], ['303', 'E', 'ILE'], ['304', 'E', 'ALA'], ['305', 'E', 'ILE'], ['306', 'E', 'TYR'], ['307', 'E', 'GLN'], ['308', 'E', 'ARG'], ['309', 'E', 'LYS'], ['310', 'E', 'VAL'], ['311', 'E', 'GLU'], ['312', 'E', 'ALA'], ['313', 'E', 'PRO'], ['314', 'E', 'PHE'], ['315', 'E', 'ILE'], ['316', 'E', 'PRO'], ['317', 'E', 'LYS'], ['318', 'E', 'PHE'], ['319', 'E', 'LYS'], ['321', 'E', 'PRO'], ['323', 'E', 'ASP'], ['324', 'E', 'THR'], ['325', 'E', 'SER'], ['326', 'E', 'ASN'], ['327', 'E', 'PHE'], ['328', 'E', 'ASP'], ['329', 'E', 'ASP'], ['330', 'E', 'TYR'], ['331', 'E', 'GLU'], ['332', 'E', 'GLU'], ['333', 'E', 'GLU'], ['334', 'E', 'GLU'], ['335', 'E', 'ILE'], ['336', 'E', 'ARG'], ['337', 'E', 'VAL'], ['338', 'E', 'SEP'], ['339', 'E', 'ILE'], ['340', 'E', 'ASN'], ['341', 'E', 'GLU'], ['342', 'E', 'LYS'], ['343', 'E', 'CYS'], ['345', 'E', 'LYS'], ['346', 'E', 'GLU'], ['347', 'E', 'PHE'], ['348', 'E', 'THR'], ['349', 'E', 'GLU'], ['350', 'E', 'PHE']]

Further, users can save the list of solvent-exposed residues to a pickle binary file, a PyMOL session or create a graph in the same fashion as above.

*************
API Reference
*************

**SERD.detect(target, surface_representation='SES', step=0.6, probe=1.4, vdw=None, ignore_backbone=True, nthreads=None, verbose=False)**

Detect solvent-exposed residues of a target biomolecule.

:Parameters:      

  * **target** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_, `pathlib.Path <https://docs.python.org/3/library/pathlib.html#pathlib.Path>`_]) – A path to PDB or XYZ file of a target biomolecular structure.

  * **surface_representation** (`Literal <https://docs.python.org/3/library/typing.html#typing.Literal>`_\["VDW", "SES", "SAS"], optional) – Surface representation. Keywords options are VDW (van der Waals), SES (Solvent Excluded Surface) or SAS (Solvent Accessible Surface), by default “SES”.

  * **step** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`float <https://docs.python.org/3/library/functions.html#float>`_, `int <https://docs.python.org/3/library/functions.html#int>`_], *optional*) – Grid spacing (A), by default 0.6.

  * **probe** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`float <https://docs.python.org/3/library/functions.html#float>`_, `int <https://docs.python.org/3/library/functions.html#int>`_], *optional*) – Probe size (A) to define SES and SAS representations, by default 1.4.

  * **vdw** (`Optional <https://docs.python.org/3/library/typing.html#typing.Optional>`_\[`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_, `pathlib.Path <https://docs.python.org/3/library/pathlib.html#pathlib.Path>`_]], *optional*) – A path to a van der Waals radii file, by default None. If None, apply the built-in van der
    Waals radii file: *vdw.dat*.

  * **ignore_backbone** (`bool <https://docs.python.org/3/library/functions.html#bool>`_, *optional*) – Whether to ignore backbone atoms (C, CA, N, O) when defining interface residues, by default True.

  * **nthreads** (`Optional <https://docs.python.org/3/library/typing.html#typing.Optional>`_\[`int <https://docs.python.org/3/library/functions.html#int>`_], *optional*) – Number of threads, by default None. If None, the number of threads is *os.cpu_count() - 1*.

  * **verbose** (`bool <https://docs.python.org/3/library/functions.html#bool>`_, *optional*) – Print extra information to standard output, by default False.

:Returns:         
  **residues** – A list of solvent-exposed residues.

:Return type:     
  `List <https://docs.python.org/3/library/typing.html#typing.List>`_\[`List <https://docs.python.org/3/library/typing.html#typing.List>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_]]

:Raises:          
  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *target* must be a string or a pathlib.Path.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *surface_representation* must be a *VDW*, *SES* or *SAS*.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *step* must be a positive real number.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *step* must be a positive real number.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *probe* must be a non-negative real number.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *probe* must be a non-negative real number.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *vdw* must be a string or a pathlib.Path.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *nthreads* must be a positive integer.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *nthreads* must be a positive integer.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *verbose* must be a boolean.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *probe* must be a positive real number, when SES or SAS is set.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *target* must be .pdb or .xyz.

.. note:: 
  
  The van der Waals radii file defines the radius values for each
  atom by residue and when not defined, it uses a generic value
  based on the atom type (see `pyKVFinder <https://github.com/LBC-LNBio/pyKVFinder>`_ package).

**SERD.read_vdw(fn=None)**

Reads van der Waals radii from .dat file.

:Parameters:      

  **fn** (`Optional <https://docs.python.org/3/library/typing.html#typing.Optional>`_\[`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_, `pathlib.Path <https://docs.python.org/3/library/pathlib.html#pathlib.Path>`_]], *optional*) – A path to a van der Waals radii file, by default None. If None, apply the built-in van der
  Waals radii file: *vdw.dat*.

:Returns:         
  **vdw** – A dictionary containing radii values.

:Return type:     
  `Dict <https://docs.python.org/3/library/typing.html#typing.Dict>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_, `Dict <https://docs.python.org/3/library/typing.html#typing.Dict>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_, `float <https://docs.python.org/3/library/functions.html#float>`_]]

:Raises:          
  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *fn* must be a string or a pathlib.Path.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – A line in *vdw* has incorrect format. The values must be double tab-separated.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – A line in *vdw* has an incorrect radius type for an atom.

.. note:: 
  
  The van der Waals radii file defines the radius values for each
  atom by residue and when not defined, it uses a generic value
  based on the atom type (see *van der Waals file template* in 
  `pyKVFinder <https://github.com/LBC-LNBio/pyKVFinder>`_ package).
  The package contains a built-in van der Waals radii file: *vdw.dat*.

**SERD.read_pdb(fn, vdw=None)**

Reads PDB file into numpy.ndarrays.

:Parameters:

  * **fn** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_, `pathlib.Path <https://docs.python.org/3/library/pathlib.html#pathlib.Path>`_]) – A path to PDB file.

  * **vdw** (`Dict <https://docs.python.org/3/library/typing.html#typing.Dict>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_, `Dict <https://docs.python.org/3/library/typing.html#typing.Dict>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_, `float <https://docs.python.org/3/library/functions.html#float>`_]], *optional*) – A dictionary containing radii values, by default None. If None, use output of *SERD.read_vdw()*.

:Returns:         
  **atomic** – A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
  and radius) for each atom.

:Return type:     
  numpy.ndarray

:Raises:          
  `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *fn* must be a string or a pathlib.Path.

.. note:: 
  
  The van der Waals radii file defines the radius values for each atom
  by residue and when not defined, it uses a generic value based on the
  atom type. The function by default loads the built-in van der Waals radii
  file: *vdw.dat*.

**SERD.get_vertices(atomic, probe=1.4, step=0.6)**

Gets 3D grid vertices.

:Parameters:

  * **atomic** (numpy.ndarray) – A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
    and radius) for each atom.

  * **probe** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`float <https://docs.python.org/3/library/functions.html#float>`_, `int <https://docs.python.org/3/library/functions.html#int>`_], *optional*) – Probe size (A), by default 4.0.

  * **step** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`float <https://docs.python.org/3/library/functions.html#float>`_, `int <https://docs.python.org/3/library/functions.html#int>`_], *optional*) – Grid spacing (A), by default 0.6.

:Returns:         
  **vertices** – A numpy.ndarray with xyz vertices coordinates
  (origin, X-axis, Y-axis, Z-axis).

:Return type:     
  numpy.ndarray

**SERD.surface(atomic, surface_representation='SES', step=0.6, probe=1.4, nthreads=None, verbose=False)**

Defines the solvent-exposed surface of a target biomolecule.

:Parameters:   

  * **atomic** (numpy.ndarray) – A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
    and radius) for each atom.

  * **surface_representation** (`Literal <https://docs.python.org/3/library/typing.html#typing.Literal>`_\["VDW", "SES", "SAS"], *optional*) – Surface representation. Keywords options are VDW (van der Waals), SES (Solvent Excluded Surface) or SAS (Solvent Accessible Surface), by default “SES”.

  * **step** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`float <https://docs.python.org/3/library/functions.html#float>`_, `int <https://docs.python.org/3/library/functions.html#int>`_], *optional*) – Grid spacing (A), by default 0.6.

  * **probe** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`float <https://docs.python.org/3/library/functions.html#float>`_, `int <https://docs.python.org/3/library/functions.html#int>`_], *optional*) – Probe size (A) to define SES and SAS representations, by default 1.4.

  * **nthreads** (`Optional <https://docs.python.org/3/library/typing.html#typing.Optional>`_\[`int <https://docs.python.org/3/library/functions.html#int>`_], *optional*) – Number of threads, by default None. If None, the number of threads is *os.cpu_count() - 1*.

  * **verbose** (`bool <https://docs.python.org/3/library/functions.html#bool>`_, *optional*) – Print extra information to standard output, by default False.

:Returns:         
  **surface** – Surface points in the 3D grid (surface[nx, ny, nz]).
  Surface array has integer labels in each positions, that are:

  * -1: solvent points;

  * 0: biomolecule points;

  * 1: solvent-exposed surface points.

  Enclosed regions are considered biomolecule points.

:Return type:     
  numpy.ndarray

:Raises:          
  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *atomic* must be a numpy.ndarray.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *atomic* has incorrect shape. It must be (n, 8).

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *surface_representation* must be a *VDW*, *SES* or *SAS*.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *step* must be a positive real number.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *step* must be a positive real number.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *probe* must be a non-negative real number.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *probe* must be a non-negative real number.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *nthreads* must be a positive integer.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *nthreads* must be a positive integer.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *verbose* must be a boolean.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *probe* must be a positive real number, when SES or SAS is set.

**SERD.interface(surface, atomic, ignore_backbone=True, step=0.6, probe=1.4, nthreads=None, verbose=False)**

Identify solvent-exposed residues based on a target solvent-exposed surface
and atomic information of a biomolecule (residue number, chain identifier, residue
name, xyz coordinates, radius).

:Parameters:      

  * **surface** (numpy.ndarray) – 

    Surface points in the 3D grid (surface[nx, ny, nz]).
    Surface array has integer labels in each positions, that are:

    * -1: solvent points;

    * 0: biomolecule points;

    * 1: solvent-exposed surface points.

    Enclosed regions are considered biomolecule points.

  * **atomic** (numpy.ndarray) – A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
    and radius) for each atom.

  * **ignore_backbone** (`bool <https://docs.python.org/3/library/functions.html#bool>`_, *optional*) – Whether to ignore backbone atoms (C, CA, N, O) when defining interface residues, by default True.

  * **step** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`float <https://docs.python.org/3/library/functions.html#float>`_, `int <https://docs.python.org/3/library/functions.html#int>`_], *optional*) – Grid spacing (A), by default 0.6.

  * **probe** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`float <https://docs.python.org/3/library/functions.html#float>`_, `int <https://docs.python.org/3/library/functions.html#int>`_], *optional*) – Probe size (A) to define SES and SAS representations, by default 1.4.

  * **nthreads** (`Optional <https://docs.python.org/3/library/typing.html#typing.Optional>`_\[`int <https://docs.python.org/3/library/functions.html#int>`_], *optional*) – Number of threads, by default None. If None, the number of threads is *os.cpu_count() - 1*.

  * **verbose** (`bool <https://docs.python.org/3/library/functions.html#bool>`_, *optional*) – Print extra information to standard output, by default False.

:Returns:         
  **residues** – A list of solvent-exposed residues.

:Return type:     
  `List <https://docs.python.org/3/library/typing.html#typing.List>`_\[`List <https://docs.python.org/3/library/typing.html#typing.List>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_]]

:Raises:          
  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *surface* must be a numpy.ndarray.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *surface* has the incorrect shape. It must be (nx, ny, nz).

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *atomic* must be a numpy.ndarray.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *atomic* has incorrect shape. It must be (n, 8).

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *ignore_backbone* must be a boolean.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *step* must be a positive real number.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *probe* must be a non-negative real number.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *probe* must be a non-negative real number.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *nthreads* must be a positive integer.

  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *nthreads* must be a positive integer.

  * `TypeError <https://docs.python.org/3/library/exceptions.html#TypeError>`_ – *verbose* must be a boolean.

**SERD.save(residues, fn='residues.pickle')**

Save list of solvent-exposed residues to binary pickle file.

:Parameters:      
  * **residues** (`List <https://docs.python.org/3/library/typing.html#typing.List>`_\[`List <https://docs.python.org/3/library/typing.html#typing.List>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_]]) – A list of solvent-exposed residues.

  * **fn** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_, `pathlib.Path <https://docs.python.org/3/library/pathlib.html#pathlib.Path>`_], *optional*) – A path to pickle file, by default “residues.pickle”

**SERD.save_session(target, residues, fn='residues.pse')**

Save a PyMOL session with the solvent-exposed residues (shown as red sticks) and
the target biomolecular structure (shown as sticks).

:Parameters:     

  * **target** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_, `pathlib.Path <https://docs.python.org/3/library/pathlib.html#pathlib.Path>`_]) – A path to PDB or XYZ file of a target biomolecular structure.

  * **residues** (`List <https://docs.python.org/3/library/typing.html#typing.List>`_\[`List <https://docs.python.org/3/library/typing.html#typing.List>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_]]) – A list of solvent-exposed residues.

  * **fn** (`Union <https://docs.python.org/3/library/typing.html#typing.Union>`_\[`str <https://docs.python.org/3/library/stdtypes.html#str>`_, `pathlib.Path <https://docs.python.org/3/library/pathlib.html#pathlib.Path>`_], *optional*) – A path to a PyMOL session file, by default “residues.pse”

**SERD.r2g(residues, atomic, selection="CB", cutoff=None, intraresidual=False)**

Create a graph from a list of solvent-exposed residues.

:Parameters:      
  
  * **residues** – A list of solvent-exposed residues.

  * **atomic** (numpy.ndarray) – A numpy array with atomic data (residue number, chain, residue name, atom name, xyz coordinates
    and radius) for each atom.

  * **selection** (`Literal <https://docs.python.org/3/library/typing.html#typing.Literal>`_\["CA", "CB", "all"], *optional*)

    Atomic selection, by default "CB". Keywords options are:
    
    * 'CA': Select alfa-carbon;

    * 'CB': Select beta-carbon, except for glycine which selects the alfa-carbon;
    
    * 'all': Select all atoms, distance between residues are the smallest distance between the atoms of these residues.

  * **cutoff** (`Optional <https://docs.python.org/3/library/typing.html#typing.Optional>`_\[`float <https://docs.python.org/3/library/functions.html#float>`_], *optional*) – A limit of distance to define an edge between two solvent-exposed residues, by default None. If None, cutoff depends on selection argument. If "CA", cutoff is 10.0. If "CB", cutoff is 8.0. If "all", cutoff is 5.0.
  
  * **intraresidual** (`bool <https://docs.python.org/3/library/functions.html#bool>`_, *optional*) – Whether to consider intraresidual contacts to create adjacency matrix, by default False.

:Returns:         
  **G** – A graph of solvent-exposed residues with edges defined by a distance smaller than the cutoff.

:Return type:     
  `networkx.classes.graph.Graph <https://networkx.org/documentation/stable/reference/classes/index.html?highlight=networkx%20classes%20graph%20graph>`_

:Raises:          
  * `ValueError <https://docs.python.org/3/library/exceptions.html#ValueError>`_ – *selection* must be *CA*, *CB*, or *all*.

.. note::
  
  Cutoff for beta-carbon is based on CAPRI round 28. For more details, refer to https://www.ebi.ac.uk/msd-srv/capri/round28/round28.html.

*******
License
*******

The software is licensed under the terms of the GNU General Public License version 3 (GPL3) and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
