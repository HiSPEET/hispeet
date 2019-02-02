
    # load environment, e.g.
    module load compiler/intel

    make distclean
    mkdir build

    # all things together
    cd build
    cmake .. 
    make -j
    
    # optional testing
    ctest

    # step by step and more informative
    cd build
    rm -rf *
    cmake ..
    make VERBOSE=1
    ctest --verbose

### Debugging

    cd build
    cmake -DCMAKE_BUILD_TYPE=Debug ..