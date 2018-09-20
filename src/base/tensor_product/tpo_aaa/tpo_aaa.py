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

operator  = 'TPO_AAA'
procedure = 'procedure(TPO_AAA_Proc)'

module    = operator.lower()
dest_proc = dest_dir + module + '__var.F'

# OpenACC device parameters
max_vec_length = 1024

#-----------------------------------------------------------------------------
# generic procedures

shutil.copy(tmpl_dir + module + '__gen.f', dest_proc)

#-----------------------------------------------------------------------------
# parametrized procedures

# read parameters ............................................................

config = dest_dir + module + '.var'

# read parameters
with open(config, 'r') as f:
    reader = csv.reader(f)
    proc_par = [ strip(p) for p in list(reader) if no_comment(p) ]

# parameters
na1 = [ entry[0] for entry in proc_par ]
na2 = [ entry[1] for entry in proc_par ]
op1 = [ entry[2] for entry in proc_par ]
op2 = [ entry[3] for entry in proc_par ]
op3 = [ entry[4] for entry in proc_par ]

# create procedures ..........................................................

tmpl_proc = tmpl_dir + module + '__par.Ft'
name_proc = []

for i in range(len(proc_par)):

    tag = na1[i] + 'x' + na2[i]

    subop_1 = tmpl_dir + 'subop_1__' + op1[i] + '.f'
    subop_2 = tmpl_dir + 'subop_2__' + op2[i] + '.f'
    subop_3 = tmpl_dir + 'subop_3__' + op3[i] + '.f'

    # derived parameters
    na1_t2 = str( int(na1[i]) // 2 * 2 )
    na1_t4 = str( int(na1[i]) // 4 * 4 )
    na1_t8 = str( int(na1[i]) // 8 * 8 )
    na2_t2 = str( int(na2[i]) // 2 * 2 )
    na2_t4 = str( int(na2[i]) // 4 * 4 )

    # OpenACC
    if min(int(na1[i]),int(na2[i])) < 8:
        vec_length  = 128
    else:
        vec_length  = 256

    num_workers = str( max_vec_length // vec_length )
    vec_length  = str( vec_length )

    # expand parametrized template
    with open(dest_proc, 'a') as f:
        f.write('\n')
        for line in fileinput.input(tmpl_proc):
            line = line.replace( '<tag>'         , tag         )
            line = line.replace( '<na1>'         , na1[i]      )
            line = line.replace( '<na1_t2>'      , na1_t2      )
            line = line.replace( '<na1_t4>'      , na1_t4      )
            line = line.replace( '<na1_t8>'      , na1_t8      )
            line = line.replace( '<na2>'         , na2[i]      )
            line = line.replace( '<na2_t2>'      , na2_t2      )
            line = line.replace( '<na2_t4>'      , na2_t4      )
            line = line.replace( '<num_workers>' , num_workers )
            line = line.replace( '<vec_length>'  , vec_length  )
            line = line.replace( '<SubOp_1>'     , subop_1     )
            line = line.replace( '<SubOp_2>'     , subop_2     )
            line = line.replace( '<SubOp_3>'     , subop_3     )
            f.write(line)

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
                    f.write('na1 == ' + na1[i] + ' .and. ')
                    f.write('na2 == ' + na2[i])
                    f.write(') ')
                    f.write('then\n')
                    f.write('    Proc => ' + name_proc[i] + '\n')
                    IF = '  else if '
                f.write('  end if\n')
