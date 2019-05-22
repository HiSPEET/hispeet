## Using FORD for documentation

### Where to find FORD

FORD resides on GitHub. For documentation see the [Wiki]
(https://github.com/Fortran-FOSS-Programmers/ford/wiki)

### How to install FORD with no root access

- First install your own pip
 
        # in some temporary directory   
        wget https://bootstrap.pypa.io/get-pip.py
        python3 get-pip.py --user

- Add the pip binary path to your PATH (best in shell start-up file)

        export PATH=~/.local/bin:$PATH

- Install FORD in your directory

        pip install ford --user

### How to use FORD for generating the HiSPEET documentation
     
     # im the HiSPEET root directory
     make doc
     
     
     
