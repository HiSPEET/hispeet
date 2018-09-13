
    module load compiler/intel

    mkdir build

    # all things together
    cd build; rm -rf * && cmake .. && make -j && ctest

    # step by step and more informative
    cd build
    rm -rf *
    cmake ..
    make VERBOSE=1
    ctest --verbose
