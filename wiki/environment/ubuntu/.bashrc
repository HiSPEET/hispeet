# Lmod
if [ -n $LMOD_CMD ]
then
  source /usr/share/lmod/lmod/init/bash
  module use $HOME/modulefiles
fi

# Path
export PATH=~/.local/bin:$PATH

