# Git Wiki

[toc]

For a more comprehensive description consider the Git [documentation](https://git-scm.com/docs) at `git-scm.com`.

## First-Time Git setup <a name="setup"></a>

For details see the  [Git Book](https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup)

### Identity

Setting your name and email address (replace mine by yours :smile:)

    git config --global user.name "Joerg Stiller"
    git config --global user.email joerg.stiller@tu-dresden.de

### Editor

Customize git to use your favorite editor, e.g. `gedit`:

    git config --global core.editor gedit

### Diff and merge tool

For file comparison and merging you may configure diff and merge tools with a graphical user interface, e.g., with Linux,

    git config --global diff.tool meld
    git config --global merge.tool meld

or, with macOS,

    git config --global diff.tool opendiff
    git config --global merge.tool opendiff

and for automatically launching these tools

    git config --global difftool.prompt false    
    git config --global mergetool.prompt false    

### Enable password caching

To get rid off retyping the password every time when contacting a remote repository you may enable password caching, e.g. on Linux

    git config --global credential.helper cache 

and on macOS

    git config --global credential.helper osxkeychain

## Starting a new repository

This is how the HiSPEET repository was created :smile:

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

## Cloning _HiSPEET_ from GitLab <a name="cloning"></a>

_HiSPEET_ can be downloaded from the GitLab server of TU Chemnitz via

     git clone https://gitlab.hrz.tu-chemnitz.de/hispeet/hispeet.git

This creates a clone of the git repository in the directory `hispeet`. Next initialize the external libraries which are incorporated as submodules:

      cd hispeet
      git submodule update --init

If the branch of the external submodule has changed, the entry in `.gitmodules` and the corresponding `cmake` module need to be updated. Following this run

      git submodule update --init --remote


## Working with branches <a name="branches"></a>

### Getting around

Show all branches

    git branch -a

Switch to the master branch

    git checkout master

Switch to another branch

    git checkout <branch-name>

### Creating a new branch

Create a new local branch

    git checkout -b <new_branch>

This command preserves uncommited changes and thus can be used to create a working branch  without spoiling the existing branch. To add the new branch to the remote repository it needs to be pushed as follows

    git push --set-upstream origin <new_branch>

### Comparing branches <a name="comparing"></a>

To compare the current branch against master branch, showing only the names of modified files

    git diff --name-status master

For details jomit the `--name-status` option. To compare a file or directory located at `path` against the master use

    git diff --name-status master -- <path>

Similarly, you can compare against another branch. Just replace `master` by the branch name as listed with `git branch -a`.  To compare any two branches:

    git diff --name-status <branch1>..<branch2>


## Making changes <a name="changes"></a>

### Updating 

Before starting to change your branch you should update the repository using

    git pull

This command downloads and merges the latest changes from the _origin_ located on the server. Typically this works smoothly. If merge conflicts appear, don't hesitate to contact an experienced developer.

### Committing

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

The next step is to commit the staged changes

    commit -m '<description of changes>'

Finally, the committed changes can pe pushed to the server

      git push

Note that for pushing you need developer access to the remote repository.

### Unstaging changes

To remove files from stage use reset HEAD where HEAD is the last commit of the current branch. This will unstage the file but maintain the modifications.

    git reset HEAD <file>

To revert the file back to the state it was in before the changes use:

    git checkout -- <file>


## Merging <a name="merging"></a>

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
    # or
    git merge --no-ff --no-commit origin/topic

Now inspect the changes

    git status
    ...
    Changes to be committed:
            modified:   foo.f
    ...

Compare/edit the changes using the difftool and commit

    git difftool HEAD foo.f
    git commit -m 'pull/merge completed'

Alternatively, you can proceed with the mergetool

    git mergetool

Note: 
To show all changes in detail omit the `name-status` option.


## Submodules <a name="submodules"></a>

Some external libraries are provided as git submodules.
For example, this is how `libxsmm` was included:

    git submodule add https://github.com/hfp/libxsmm.git external/libxsmm
    git commit -m 'added libxsmm submodule'

Note that the submodules must be manually intialized, see [Cloning](#cloning) above.
To update the submodules from external repositories use

    git submodule update --recursive --remote

For more details on submodules confer the [Git Book](https://git-scm.com/book/en/Git-Tools-Submodules).

