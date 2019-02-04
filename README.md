# LAMBReX :sheep::crown:
Lattice Boltzmann code built on AMReX

## Build Instructions

### AMReX
**Important**: The library that AMReX builds with GNU Make is subtly different from the one that it builds with cmake. LAMBReX presently supports only the GNU Make version, though this is expected to change to cmake in the near future.

To build the static AMReX library with GNU Make in way that is compatible with the current version of LAMBReX, please follow the commands given in BUILDNOTES.txt.

### LAMBReX
At present LAMBReX is built using GNU Make, though this will be changed to cmake in the near future. Run `make install` in the src directory to build and install the static LAMBReX library to the liblambrex directory.

### Examples
At present there is only one example code, a simple simulation of a "pulse" of enhanced density on a single plane, with fully periodic boundary conditions. This calculation was chosen for its simplicity, and the output was compared against the original code on which the Lattice Boltzmann calculation part of LAMBReX is based.

Again, the example can be built by using GNU Make. Simply run `make` in the examples directory.
