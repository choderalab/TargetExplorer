#!/usr/bin/env python
#
# Initialize project directory by creating necessary subdirectories and a project metadata file.
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 11 Mar 2014

import os, argparse
import TargetExplorer

# ========
# Parse command-line arguments
# ========

argparser = argparse.ArgumentParser(description='Initialize a TargetExplorer database project by creating the necessary subdirectories.')
argparser.add_argument('--project_dir', type=str, default='.', help='(Default: ".") Optionally provide a directory path in which to initialize the database project.')
args = argparser.parse_args()

project_dir = os.path.abspath(args.project_dir)

TargetExplorer.core.init_project(project_dir)

