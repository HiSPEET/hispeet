# Git

Online ressources

  * [Documentation](https://git-scm.com/docs) at `git-scm.com`

## First-Time Git Setup

For details see [here](https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup)

### Identity

    git config --global user.name "Joerg Stiller"
    git config --global user.email joerg.stiller@tu-dresden.de

### Editor

    # Change from default editor to your own choice, e.g. 
    git config --global core.editor bbedit

### Diff and merge tool

    # Command line
    git config --global diff.tool vimdiff
    git config --global merge.tool vimdiff
    
    # MacOS - GUI
    git config --global diff.tool opendiff
    git config --global merge.tool opendiff
    
    # Linux - GUI alternatives: gvim, meld, kdiff3, e.g.
    git config --global diff.tool meld
    git config --global merge.tool meld

Do not prompt before launching the diff and merge tools

    git config --global difftool.prompt false    
    git config --global mergetool.prompt false    

### Enable password caching

To get rid off retyping the password every time when contacting a remote repository you may enable password caching. Note, however, that cloning the repository using the HTTPS protocol often provides a more convenient approach. 

With _MacOS_ use the osxkeychain credential helper 

    git config --global credential.helper osxkeychain

Enable password chaching (900 seconds are the default;)

    git config credential.helper cache
    # or, using a longer caching period
    git config credential.helper 'cache --timeout=1200'

### How to start a new Git repository

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

     git clone https://scm.fusionforge.zih.tu-dresden.de/anonscm/git/hispeet/hispeet.git

Alternatively, you may pick your version of the `clone` command from 
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

### Create a new branch

__Note:__ This may fail when using the HTTPS protocol :(

Create a new (local) branch

    git checkout -b <new_branch>

Push the branch on the remote repository

    git push origin <new_branch>

### Making changes and update your branch

Unless told otherwise use the `master` branch.
But you can change to another branch using

    git checkout <branch-name>

Unless you made changes you can update your local repo using

    make distclean
    git pull

in the `HiSPEET` directory.
Before pulling, all changes to the current branch should be committed or reverted.

### Staging and committing

To save your work, visit the changes using

    git status

To stage new or modified file or a new directory

    git add <file>
    git add <dir>

Of course, files and directories can also be removed

    git rm <file>
    git rm -r <dir>

To remove a file from the repository but keep it on disk, use `--cache`

    git rm <file> --cache

Finally, commit the staged changes

    commit -m '<description of changes>'

### Unstaging

To remove files from stage use reset HEAD where HEAD is the last commit of the current branch. This will unstage the file but maintain the modifications.

    git reset HEAD <file>

To revert the file back to the state it was in before the changes use:

    git checkout -- <file>

### Comparing branches

To compare the current branch against master branch, showing only the names of modified files

    git diff --name-status master

For details jomit the `--name-status` option. To compare a file or directory located at `path` against the master use

    git diff --name-status master -- <path>

Similarly, you can compare against another branch. Just replace `master` by the branch name as listed with `git branch -a`.  To compare any two branches:

    git diff --name-status <branch1>..<branch2>


## Merging

### Basic merging

Merge branch `topic` into current branch

    git merge topic

Auto-resolve conflicts

    # favor the current branch
    git merge topic -X ours
    
    # favor other branch (topic)
    git merge topic -X theirs

### Advanced merging

Disable fast-forward merging and do not create a merge commit.
This is useful if auto-merging cannot be trusted.
Works with `pull` as well as with `merge`

    git pull --no-ff --no-commit

Now inspect the changes

    git status
    ...
    Changes to be committed:
            modified:   foo.f
    ...

Compare/edit the changes using the difftool and commit

    git difftool HEAD foo.f
    git commit -m 'pull completed'


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

