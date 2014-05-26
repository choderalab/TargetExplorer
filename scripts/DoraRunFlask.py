#!/usr/bin/env python
import os, sys
this_dir = os.path.dirname(os.path.realpath(__file__))
target_explorer_basedir = os.path.realpath(os.path.join(this_dir, '..'))
sys.path.append(target_explorer_basedir)

import app_master
app_master.app.run(debug=True)
