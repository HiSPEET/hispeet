#! /usr/bin/env python3

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

operator = 'CART__TPO_Div'
module = operator.lower()

procedure = 'procedure(TPO_Div_Proc)'

# OpenACC device parameters
max_vec_length = 1024

#-----------------------------------------------------------------------------
# generic procedure

shutil.copy(tmpl_dir + module + '__gen.f', dest_dir)

#-----------------------------------------------------------------------------
# parametrized procedures

# read parameters ............................................................

config = dest_dir + module + '.var'

# read parameters
with open(config, 'r') as f:
    reader = csv.reader(f)
    proc_par = [ strip(p) for p in list(reader) if no_comment(p) ]

# shorthands
np  = [ entry[0] for entry in proc_par ]
op1 = [ entry[1] for entry in proc_par ]
op2 = [ entry[2] for entry in proc_par ]
op3 = [ entry[3] for entry in proc_par ]

# create procedures ..........................................................

tmpl_proc = tmpl_dir + module + '__par.Ft'
dest_proc = dest_dir + module + '__par.F'
name_proc = []

for i in range(len(proc_par)):

    tag = np[i]

    subop_1i = shared_dir + 'subop_1i__' + op1[i] + '.f'
    subop_2a = shared_dir + 'subop_2a__' + op2[i] + '.f'
    subop_3a = shared_dir + 'subop_3a__' + op3[i] + '.f'

    # derived parameters
    np_t2 = str( int(np[i]) // 2 * 2 )
    np_t4 = str( int(np[i]) // 4 * 4 )
    np_t8 = str( int(np[i]) // 8 * 8 )

    # OpenACC
    if int(np[i]) < 8:
        vec_length  = 128
    else:
        vec_length  = 256

    num_workers = str( max_vec_length // vec_length )
    vec_length  = str( vec_length )

    # expand parametrized template
    with open(dest_proc, 'a') as f:
        for line in fileinput.input(tmpl_proc):
            line = line.replace( '<tag>'         , tag         )
            line = line.replace( '<np>'          , np[i]       )
            line = line.replace( '<np_t2>'       , np_t2       )
            line = line.replace( '<np_t4>'       , np_t4       )
            line = line.replace( '<np_t8>'       , np_t8       )
            line = line.replace( '<num_workers>' , num_workers )
            line = line.replace( '<vec_length>'  , vec_length  )
            line = line.replace( '<SubOp_1i>'    , subop_1i    )
            line = line.replace( '<SubOp_2a>'    , subop_2a    )
            line = line.replace( '<SubOp_3a>'    , subop_3a    )
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
                f.write('\n  select case(np)\n')
                for i in range(len(proc_par)):
                    f.write('  case(' + np[i] + ')\n')
                    f.write('    Proc => ' + name_proc[i] + '\n')
                f.write('  end select\n')
