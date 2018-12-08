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

operator  = 'CART__TPO_Elliptic_VI'
procedure = 'procedure(TPO_Elliptic_VI_Proc)'

module    = operator.lower()
dest_proc = dest_dir + module + '__var.F'

# OpenACC device parameters
max_vec_length = 1024

#-----------------------------------------------------------------------------
# generic procedures

gen_proc = [ tmpl_dir + module + '__gen.f' ] #,
#             tmpl_dir + module + '__gen_2.f' ,
#             tmpl_dir + module + '__gen_3.f' ,
#             tmpl_dir + module + '__gen_4.f' ]

with open(dest_proc, 'a') as f:
    for proc in gen_proc:
        with open(proc, 'r') as p: f.write(p.read())

#x# #-----------------------------------------------------------------------------
#x# # parametrized procedures
#x#
#x# # read parameters ............................................................
#x#
#x# config = dest_dir + module + '.var'
#x#
#x# # read parameters
#x# with open(config, 'r') as f:
#x#     reader = csv.reader(f)
#x#     proc_par = [ strip(p) for p in list(reader) if no_comment(p) ]
#x#
#x# # shorthands
#x# np  = [ entry[0] for entry in proc_par ]
#x# op1 = [ entry[1] for entry in proc_par ]
#x# op2 = [ entry[2] for entry in proc_par ]
#x# op3 = [ entry[3] for entry in proc_par ]
#x#
#x# # create procedures ..........................................................
#x#
#x# tmpl_proc = tmpl_dir + module + '__par.Ft'
#x# name_proc = []
#x#
#x# for i in range(len(proc_par)):
#x#
#x#     tag = np[i]
#x#
#x#     subop_1a = shared_dir + 'subop_1a__' + op1[i] + '.f'
#x#     subop_2a = shared_dir + 'subop_2a__' + op2[i] + '.f'
#x#     subop_3a = shared_dir + 'subop_3a__' + op3[i] + '.f'
#x#
#x#     # derived parameters
#x#     np_t2 = str( int(np[i]) // 2 * 2 )
#x#     np_t4 = str( int(np[i]) // 4 * 4 )
#x#     np_t8 = str( int(np[i]) // 8 * 8 )
#x#
#x#     # OpenACC
#x#     if int(np[i]) < 8:
#x#         vec_length  = 128
#x#     else:
#x#         vec_length  = 256
#x#
#x#     num_workers = str( max_vec_length // vec_length )
#x#     vec_length  = str( vec_length )
#x#
#x#     # expand parametrized template
#x#     with open(dest_proc, 'a') as f:
#x#         for line in fileinput.input(tmpl_proc):
#x#             line = line.replace( '<tag>'         , tag         )
#x#             line = line.replace( '<np>'          , np[i]       )
#x#             line = line.replace( '<np_t2>'       , np_t2       )
#x#             line = line.replace( '<np_t4>'       , np_t4       )
#x#             line = line.replace( '<np_t8>'       , np_t8       )
#x#             line = line.replace( '<num_workers>' , num_workers )
#x#             line = line.replace( '<vec_length>'  , vec_length  )
#x#             line = line.replace( '<SubOp_1a>'    , subop_1a    )
#x#             line = line.replace( '<SubOp_2a>'    , subop_2a    )
#x#             line = line.replace( '<SubOp_3a>'    , subop_3a    )
#x#             f.write(line)
#x#
#x#         f.write('\n')
#x#
#x#     name_proc.append(operator + '__' +  tag)

#-----------------------------------------------------------------------------
# operator module

tmpl_module = tmpl_dir + module + '.ft'
dest_module = dest_dir + module + '.f'

with open(dest_module, 'w') as f:
    for line in fileinput.input(tmpl_module):

        f.write(line)

#x#         if len(proc_par) > 0:
#x#
#x#             # external statements ............................................
#x#
#x#             if '! external procedures ...' in line:
#x#                 for i in range(len(proc_par)):
#x#                     f.write('\n  ' + procedure + ' :: ' + name_proc[i])
#x#
#x#             # assignments to parametrized procedures .........................
#x#
#x#             if '! parametrized procedures ...' in line:
#x#                 f.write('\n  select case(np)\n')
#x#                 for i in range(len(proc_par)):
#x#                     f.write('  case(' + np[i] + ')\n')
#x#                     f.write('    Proc => ' + name_proc[i] + '\n')
#x#                 f.write('  end select\n')
#x#