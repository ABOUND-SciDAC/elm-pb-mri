# ELM-pb

Multi-rate version of ELM-pb 3-field model.

# Building

This repository includes BOUT++ as a submodule that will be built by
default. To build elm-pb-mri with the submodule version of BOUT++,
downloading and linking with SUNDIALS:

    $ cmake . -B build -DBOUT_DOWNLOAD_SUNDIALS=ON
    $ cmake --build build

On MacOS the following configuration has worked:

    $ cmake . -B build -DBOUT_DOWNLOAD_SUNDIALS=ON -DBOUT_ENABLE_BACKTRACE=Off -DBUILD_SHARED_LIBS=Off -DBOUT_USE_NLS=Off -DBOUT_USE_UUID_SYSTEM_GENERATOR=Off

If you have already installed BOUT++ and want to use that rather than
configure and build BOUT++ again, set `ELMPB_BUILD_BOUT` to `OFF` and pass
CMake the path to the BOUT++ `build` directory e.g.

    $ cmake . -B build -DELMPB_BUILD_BOUT=OFF -DCMAKE_PREFIX_PATH=$HOME/BOUT-dev/build

# Running

After successfully compiling, the `build` subdirectory should contain
the `elm-pb-mri` executable.

Execute with the following:

    $ ./elm-pb-mri -d examples/data
