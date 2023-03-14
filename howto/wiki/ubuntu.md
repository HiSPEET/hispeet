## Ubuntu

Installing _HiSPEET_ on Ubuntu 22 requires the packages

- `git`
- `cmake`
- `gcc`
- `g++`
- `gfortran`
- `graphviz`
- `libopenblas-dev`
- `libfftw3-dev`
- `libmetis-dev`
- `libparmetis-dev`
- `libvtk7-dev`
- `lmod`
- `make`
- `openmpi-bin`
- `parmetis-doc`
- `python3`
- `python3-matplotlib`
- `python3-numpy`
- `paraview`

Do not install `libvtk9-dev`, which is incompatible with `paraview`

For installing a package you may use the command

    sudo apt-get install <package name>

Additionally, you may want to install the `ford` documentation package

    pip install ford


