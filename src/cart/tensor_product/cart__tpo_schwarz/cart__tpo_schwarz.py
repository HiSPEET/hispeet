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
shared_dir = args.shared + '/'

operator  = 'CART__TPO_Schwarz'
procedure = 'procedure(TPO_Schwarz_Proc)'

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
n1  = [ entry[0] for entry in proc_par ]
n2  = [ entry[1] for entry in proc_par ]
n3  = [ entry[2] for entry in proc_par ]
op1 = [ entry[3] for entry in proc_par ]
op2 = [ entry[4] for entry in proc_par ]
op3 = [ entry[5] for entry in proc_par ]

# create procedures ..........................................................

tmpl_proc = tmpl_dir + module + '__par.Ft'
name_proc = []

for i in range(len(proc_par)):

    tag = n1[i] + 'x' + n2[i] + 'x' + n3[i]

    subop_1id = shared_dir + 'subop_1id__' + op1[i] + '.f'
    subop_1i  = shared_dir + 'subop_1i__'  + op1[i] + '.f'
    subop_2i  = shared_dir + 'subop_2i__'  + op2[i] + '.f'
    subop_3i  = shared_dir + 'subop_3i__'  + op3[i] + '.f'

    # OpenACC
    if min(int(n1[i]), int(n2[i]), int(n3[i])) < 8:
        vec_length  = 128
    else:
        vec_length  = 256

    num_workers = str( max_vec_length // vec_length )
    vec_length  = str( vec_length )

    # expand parametrized template
    with open(dest_proc, 'a') as f:
        for line in fileinput.input(tmpl_proc):
            line = line.replace('<tag>'          , tag         )
            line = line.replace('<n1>'           , n1[i]       )
            line = line.replace('<n2>'           , n2[i]       )
            line = line.replace('<n3>'           , n3[i]       )
            line = line.replace( '<num_workers>' , num_workers )
            line = line.replace( '<vec_length>'  , vec_length  )
            line = line.replace( '<SubOp_1id>'   , subop_1id   )
            line = line.replace( '<SubOp_1i>'    , subop_1i    )
            line = line.replace( '<SubOp_2i>'    , subop_2i    )
            line = line.replace( '<SubOp_3i>'    , subop_3i    )
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
                    f.write('n1 == ' + n1[i] + ' .and. ')
                    f.write('n2 == ' + n2[i] + ' .and. ')
                    f.write('n3 == ' + n3[i])
                    f.write(') ')
                    f.write('then\n')
                    f.write('    Proc => ' + name_proc[i] + '\n')
                    IF = '  else if '
                f.write('  end if\n')
