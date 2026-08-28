# Editors
To configure editors for use with `HiSPEET` we recommend

1. to install the Fortran Language Server, and
2. to associate the file extensions `.f` and `.F` with free source form.

## Fortran Language Server

The [Fortran Language Server](https://pypi.org/project/fortls) enables syntax checking. It can be installed system-wide using
     pip install fortls
or only in user space with
     pip install --user fortls
Depending on your system, other installation options may exist.

## Vim
To set free source form for all Fortran files, add the following lines to your `.vimrc`

     let b:fortran_fixed_source=0
     let b:fortran_free_source=1
     let b:fortran_more_precise=1
     let b:fortran_do_enddo=1

## Visual Studio Code
First, go to menu `Code > Preferences > Extensions` and install `Modern Fortran`. To associate the extensions `.f` and `.F` with free source form, open or create `~/.vscode/settings.json` and insert the following lines

     "files.associations": {
           "*.f": "FortranFreeForm",
           "*.F": "FortranFreeForm"
        }

Optionally, you can specify compiler dependent options. For example, supposing you use the `gfortran` located in `/usr/local/bin`:

        "[FortranFreeForm]": {
    
          },
        "fortran.linter": "gfortran",
        "fortran.linter.executablePath": "/usr/local/bin/gfortran",
        "fortran.linter.includePaths": [ "${workspaceFolder}/**" ],
        "fortran.linter.extraArgs": ["-ffree-form", "-std=f2008"],
        "fortran.fortls.preprocessor.suffixes": [ ".F"]

Please note that the options alltogether are enclosed in curly brackets, i.e.

      {
          // insert options here ...
      }
