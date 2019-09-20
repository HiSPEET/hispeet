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

## Cloning <a name="cloning"></a>

Clone HiSPEET into directory hispeet (will be created if not yet existing) 

     git clone https://<uid>@scm.fusionforge.zih.tu-dresden.de/authscm/<uid>/git/hispeet/hispeet.git hispeet

where `<uid>`is your Fusionforge user ID. 
You may pick your version of the `clone` command from 
[SCM menu](https://fusionforge.zih.tu-dresden.de/scm/?group_id=747):

  *  go to "Developer Access", change to "via smart HTTP",
  *  copy the command, and
  *  add the target directory

Then, change to the repo and initialize the git-submodules, e.g.

      cd hispeet
      git submodule update --init


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


## Submodules

Some external libraries are provided as git submodules.
For example, this is how `libxsmm` was included:

    git submodule add https://github.com/hfp/libxsmm.git external/libxsmm
    git commit -m 'added libxsmm submodule'

Note that the submodules must be manually intialized, see [Cloning](#cloning) above.
Also you may wish to update the submodules from the external repository, e.g.

    git submodule update

 
