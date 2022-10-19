Requirements
============

If you don’t have Python v3 and/or SWIG yet, the installation procedure differs depending on the operating system.

Python and SWIG can be installed with conda, or with a package manager on Linux and macOS.

Package managers
----------------

On Linux:

.. code-block:: bash
    
  sudo apt install python3 swig

On macOS:

.. code-block:: bash
    
  brew install python3 swig

.. note:: 

  Users can use their preferred package manager to install SWIG and Python v3.

Conda
-----

If you use conda, you can install Python v3 and SWIG from the defaults channel:

.. code-block:: bash
    
  # Use an environment rather than install in base environement (recommended)
  conda create -n myenv python=3
  conda activate myenv
  # The SWIG install command
  conda install swig

Installation
============

The prerequisites for installing SERD is SWIG, Python v3 and PyMOL v2.

To install the latest release on `PyPI <https://pypi.org/project/SERD>`_, 
run:

::

  pip install SERD

Or to install the latest developmental version, run:

::

  git clone https://github.com/jvsguerra/SERD.git
  pip install SERD
