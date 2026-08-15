# Build Wiki

## Basics

Initialize the computing environment as described in the [environment](./environment) wiki. Then proceed as follows:

- Go to the _HiSPEET_ root directory

- Optionally, remove old build directories

      make distclean

- Create the build directory, e.g.

      mkdir build

- Generate the configuratiion and and accomplish a complete build

      cd build
      cmake .. 
      make -j

  In case of trouble try a sequential build, i.e.

      make

## Optional testing

**currently not available**

Following to build

    ctest

Step by step and more informative

    cd build
    rm -rf *
    cmake ..
    make VERBOSE=1
    ctest --verbose

## OpenMP

Parts of the code are prepared for acceleration with OpenMP.
First build the libraries, e.g.

    mkdir build-omp
    cd build-omp
    cmake -D OpenMP=1 ..
    make -j

Now try the Cartesian tensor-product operators

    cd develop/cart/tensor-product

For example, you may test the divergence operator

    # increasing number of threads from 1 to 4
    export OMP_NUM_THREADS=1; ./validate__cart__tpo_div
    export OMP_NUM_THREADS=2; ./validate__cart__tpo_div
    export OMP_NUM_THREADS=3; ./validate__cart__tpo_div
    export OMP_NUM_THREADS=4; ./validate__cart__tpo_div

For reaching a reasonable speedup, the number of elements is important. You can see this by changing the parameter `ne` in
`validate__cart__tpo_div.prm` and repeating the test.

## Debugging

The first step in tracking down arrows is to enable runtime checks. This is done by choosing the build type `Check`: 

    cmake -D CMAKE_BUILD_TYPE=Check ..
    make

 The next option is to enable debugging:

    cmake -D CMAKE_BUILD_TYPE=Debug ..
    make

Then change to the working directory with the executable and invoke the debugger, e.g.

* `gdb` if you are using `gfortran` (`ggdb` with MacPorts under MacOS)
* `gdb-aa` when using `ifort`