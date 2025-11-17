Command Line Interface
=======================

.. code-block:: text

  usage: irmsd [-h] [--cn] [--rot] [--canonical] [--rmsd] [--irmsd] [--heavy]
               [--inversion {on,off,auto}] [-o OUTPUT]
               structures [structures ...]

  CLI to read an arbitrary number of structures with ASE and run selected
  analysis commands on them.

  positional arguments:
    structures            Paths to structure files (e.g. .xyz, .pdb, .cif). You
                          can pass many.

  options:
    -h, --help            show this help message and exit
    --cn                  Calculate coordination numbers for each structure and
                          print them as numpy arrays.
    --rot                 Calculate the rotational constants.
    --canonical           Calculate the canonical identifiers.
    --rmsd                Calculate the Cartesian RMSD between two given
                          structures via a quaternion algorithm.
    --irmsd               Calculate the invariant Cartesian RMSD between two
                          given structures.
    --heavy               When calculating RMSD or canonical atom identifier,
                          consider only heavy atoms.
    --inversion {on,off,auto}
                          Control coordinate inversion in irmsd runtypes: 'on',
                          'off', or 'auto' (default: auto).
    -o, --output OUTPUT   Output file name (optional). If not provided, nothing
                          is written.
