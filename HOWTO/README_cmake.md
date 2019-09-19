
Load environment, e.g.
  
    module load hispeet/intel

Create build directory

    make distclean
    mkdir build

Optionally, set tuning flags

    #export TUNING_FLAG=simple
    #export TUNING_FLAG=intel

Build (complete)

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