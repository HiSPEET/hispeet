## Prepare your home

      # insert 'environment/taurus/.bashrc_all' into your bash startup file
      cp -r environment/taurus/modulefiles ~

## Clone HiSPEET from GitLab 

      mkdir HiSPEET
      cd HiSPEET/
      git clone https://gitlab.hrz.tu-chemnitz.de/hispeet/hispeet.git .
      git submodule update --init

## Build HiSPEET

With Intel

     mkdir build-intel
     cd build-intel
     ml hispeet/intel-impi-2021.6.0
     cmake ..
     make

or with GNU

     mkdir build-gnu
     cd build-gnu
     ml hispeet/gnu-foss-2021b.lua
     cmake ..
     make


