# LAMBReX :sheep::crown:
Lattice Boltzmann code built on AMReX

## Build Instructions

### AMReX
**Important**: The library that AMReX builds with GNU Make is subtly different from the one that it builds with cmake. LAMBReX now supports only the cmake version. Additionally AMReX is updated *regularly*, for LAMBReX is attempting to move with it, so minimum version 19.02 is `REQUIRED`.

To build the static AMReX library with cmake in way that is compatible with the current version of LAMBReX, please follow the commands given in BUILDNOTES.txt.

### LAMBReX
From the home directory create a build directory and cd to it -- for example `mkdir build` followed by `cd build`. Then simply run `cmake ../` followed by `make`. `liblambrex.a` can be found in `${LAMBReX_HOME}/liblambrex`, the binary for the example calculation can be found in `${LAMBREX_HOME}/bins`.

### Examples
At present there is only one example code, a simple simulation of a "pulse" of enhanced density on a single plane, with fully periodic boundary conditions. This calculation was chosen for its simplicity, and the output was compared against the original code on which the Lattice Boltzmann calculation part of LAMBReX is based.
