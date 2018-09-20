#!/usr/bin/env python3

import os
import sys
import argparse
import filecmp
import shutil

#-----------------------------------------------------------------------------
# Initialization

parser = argparse.ArgumentParser()
parser.add_argument('--op'  , help='operator name')
parser.add_argument('--tmpl', help='template directory')
parser.add_argument('--dest', help='destination directory')

operator = parser.parse_args().op
tmpl_dir = parser.parse_args().tmpl + '/'
dest_dir = parser.parse_args().dest + '/'
config   = operator.lower() + '.var'

if not os.path.isdir(tmpl_dir):
    print('Template directory "', tmpl_dir, '" does not exist')
    sys.exit(1)

if not os.path.isdir(dest_dir):
    print('Destination directory "', dest_dir, '" does not exist')
    sys.exit(1)

dest_config    = dest_dir + config
default_config = tmpl_dir + config + '.default'

#-----------------------------------------------------------------------------
# Identify target configuration

tuning_flag = os.getenv('TUNING_FLAG')
tuning_dir  = os.getenv('TUNING_DIR')

if tuning_flag:
    tuning_flag  = tuning_flag.strip()
    tuned_config = config + '.' + tuning_flag

    if tuning_dir:
        tuning_dir   = tuning_dir.strip() + '/'
        tuned_config = tuning_dir + tuned_config
    else:
        tuned_config = tmpl_dir   + tuned_config

    if os.path.exists(tuned_config):
        target_config = tuned_config
    else:
        target_config = default_config

else:
    target_config = default_config

#-----------------------------------------------------------------------------
# Update configuration if changed


if os.path.isfile(dest_config):

    # replace configuration if existing, but different
    if not filecmp.cmp(target_config, dest_config):
        shutil.copy(target_config, dest_config)

else:
      shutil.copy(target_config, dest_config)
