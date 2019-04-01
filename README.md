# LAMBReX :sheep::crown:
Lattice Boltzmann code built on AMReX

## Build Instructions

### LAMBReX
From the home directory create a build directory and cd to it -- for example `mkdir build` followed by `cd build`. Then simply run `conan install ..`, `cmake ..`, and `make`. `liblambrex.a` can be found in `build/lib`, the binaries for the example calculation and tests can be found in `build/bin`.

### Dependencies

#### Conan
[Conan](https://conan.io/) is a C/C++ package manager. The recommended way to install it is with `pip`, simply using `pip install conan`. If you don't have and don't want pip installed, alternatives can be found [here](https://docs.conan.io/en/latest/installation.html).

#### AMReX
AMReX is a C++ adaptive mesh refinement library. The source can be downloaded from the [AMReX GitHub](https://github.com/AMReX-Codes/amrex).

**Important**: The library that AMReX builds with GNU Make is subtly different from the one that it builds with cmake. LAMBReX now supports only the cmake version. Additionally AMReX is updated *regularly*. LAMBReX is attempting to move with it, so minimum version 19.03 is `REQUIRED`. Additionally version 19.03+ of AMReX depends on [Blitz++](https://github.com/blitzpp/blitz) being installed.

To build the static AMReX library with cmake in way that is compatible with the current version of LAMBReX, please follow the commands given in (or run) [amrex_cmake.sh](https://github.com/otbrown/LAMBReX/blob/master/amrex_cmake.sh).

#### Catch2
The C++ Automated Tests in a Header (Catch2) test framework is used for testing LAMBReX. It can be found at [Catch2](https://github.com/catchorg/Catch2). For compatibility with LAMBReX, system-wide installation using cmake is recommended, as described [here](https://github.com/catchorg/Catch2/blob/master/docs/cmake-integration.md#installing-catch2-from-git-repository).

#### Examples
At present there is only one example code, a simple simulation of a "pulse" of enhanced density on a single plane, with fully periodic boundary conditions. This calculation was chosen for its simplicity, and the output was compared against the original code on which the Lattice Boltzmann calculation part of LAMBReX is based.
