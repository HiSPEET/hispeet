#!/usr/bin/env python3

import os
import pathlib
import argparse

#-----------------------------------------------------------------------------
#> Returns True if the given line record does not present a comment

def no_comment(line):
    record = ''.join(line).strip()
    if (len(record) > 0):
        return record[0].isalnum()
    else:
        return False

#-----------------------------------------------------------------------------
#> Strips every string in the given list

def strip(strings):
    return [ s.strip() for s in strings ]


#-----------------------------------------------------------------------------
#> Strips every string in the given list

def new_empty_file(name):
    if os.path.exists(name):
      os.remove(name)
    pathlib.Path(name).touch()


#-----------------------------------------------------------------------------
#> Get template, destination and shared routines directories from sys.argv
#>
#> Example script 'tpo.py'
#>
#>      #! /usr/bin/env python3
#>      import sys
#>      from tpo_utilities import *
#>      args   = get_args(sys.argv)
#>      tmpl   = args.tmpl
#>      dest   = args.dest
#>      shared = args.shared
#>      ...
#>
#> to be called as
#>
#>     ./tpo.py --tmpl ${TMPL} --dest ${DEST}
#>
#> Non-required arguments can be omitted.

def get_args(argv=None):
  parser = argparse.ArgumentParser()
  parser.add_argument('--tmpl'  , help='template directory')
  parser.add_argument('--dest'  , help='destination directory')
  parser.add_argument('--shared', help='shared routines directory')
  return parser.parse_args(argv[1:])
