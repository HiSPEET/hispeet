## Ubuntu configuration

- Insert the contents of `.bashrc` into you bash startup script

- Copy the folder `modulefiles` to your own module directory or,
  if you still do not have one, to your home directory, i.e.

      cp -r modulefiles ~ 
  
- When you open a new terminal you can initialize the `HiSPEET` environment using
  
      module load hispeet/gnu
  
- Now you are ready to start `cmake` from your build directory

