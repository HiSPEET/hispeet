## How it started

     mkdir HiSPEET
     cd HiSPEET

     git init
     mkdir doc
     mkdir src
     ..
     git add doc src
     git commit

     # stage the modified and deleted files
     git add -u

     # add and commit only the modified and deleted files.
     git commit -a

## Clone a local git repository from Fusionforge

TBD

## Branches

### Getting around

Show all branches

    git branch -a

Switch to your branch

    git checkout <branch-name>


### Making changes and update your branch

Unless told otherwise use the `master` branch.
But you can change to another branch using

    git checkout <branch-name>

Unless you made changes you can update your local repo using

    make distclean
    git pull

in the `HiSPEET` directory.


### Comparing branches

To compare the current branch against master branch, showing only the names of modified files

    git diff --name-status master

To compare any two branches:

$ git diff --name-status <branch1>..<branch2>

Note: 
To show all changes in detail omit the `name-status` option.
