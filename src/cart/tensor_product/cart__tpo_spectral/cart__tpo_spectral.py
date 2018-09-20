#!/usr/bin/env python3

import os
import sys
import shutil
import fileinput
import csv

# own modules
from tpo_utilities import *

#-----------------------------------------------------------------------------
# initialization

args = get_args(sys.argv)
tmpl_dir = args.tmpl + '/'
dest_dir = args.dest + '/'
#shared_dir = args.shared + '/'

operator  = 'CART__TPO_Spectral'
procedure = 'procedure(TPO_Spectral_Proc)'

module    = operator.lower()
dest_proc = dest_dir + module + '__var.F'

# OpenACC device parameters
max_vec_length = 1024

#-----------------------------------------------------------------------------
# generic procedure

shutil.copy(tmpl_dir + module + '__gen.f', dest_proc)

#-----------------------------------------------------------------------------
# parametrized procedures

# read parameters ............................................................

config = dest_dir + module + '.var'

# read parameters
with open(config, 'r') as f:
    reader = csv.reader(f)
    proc_par = [ strip(p) for p in list(reader) if no_comment(p) ]

# shorthands
na  = [ entry[0] for entry in proc_par ]
nb  = [ entry[1] for entry in proc_par ]
nc  = [ entry[2] for entry in proc_par ]
var = [ entry[3] for entry in proc_par ]

# create procedures ..........................................................

tmpl_proc = tmpl_dir + module + '__par.Ft'
name_proc = []

for i in range(len(proc_par)):

    tag = na[i] + 'x' + nb[i] + 'x' + nc[i]

    # expand parametrized template
    with open(dest_proc, 'a') as f:
        for line in fileinput.input(tmpl_proc):
            line = line.replace('<tag>', tag)
            line = line.replace('<na>', na[i])
            line = line.replace('<nb>', nb[i])
            line = line.replace('<nc>', nc[i])
            f.write(line)

        f.write('\n')

    name_proc.append(operator + '__' +  tag)

#-----------------------------------------------------------------------------
# operator module

tmpl_module = tmpl_dir + module + '.ft'
dest_module = dest_dir + module + '.f'

with open(dest_module, 'w') as f:
    for line in fileinput.input(tmpl_module):

        f.write(line)

        if len(proc_par) > 0:

            # external statements ............................................

            if '! external procedures ...' in line:
                for i in range(len(proc_par)):
                    f.write('\n  ' + procedure + ' :: ' + name_proc[i])

            # assignments to parametrized procedures .........................

            if '! parametrized procedures ...' in line:
                IF = '\n  if '
                for i in range(len(proc_par)):
                    f.write(IF)
                    f.write('(')
                    f.write('na == ' + na[i] + ' .and. ')
                    f.write('nb == ' + nb[i] + ' .and. ')
                    f.write('nc == ' + nc[i])
                    f.write(') ')
                    f.write('then\n')
                    f.write('    Proc => ' + name_proc[i] + '\n')
                    IF = '  else if '
                f.write('  end if\n')
