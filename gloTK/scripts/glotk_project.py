#!/usr/bin/env python
# gloTK - Genomes of Luminous Organisms Toolkit
# Copyright (c) 2015-2016 Darrin Schultz. All rights reserved.
#
# This file is part of gloTK.
#
# GloTK is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GloTK is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GloTK.  If not, see <http://www.gnu.org/licenses/>.

"""

title: glotk-project
authr: darrin schultz

This program:
1. Reads in a meraculous config file and initializes a glotk assembly project in
   the current directory.

Usage:
glotk-project --inputConfig <location of Meraculous config> --genus pleu --species bach
OR
gglotk-project -i <location of Meraculous config> -g pleu -s bach
"""

#import things for rest of program
import argparse
import copy
import os
import shutil
import sys


#import gloTK stuff
from gloTK import ConfigParse
import gloTK.utils

#This class is used in argparse to expand the ~. This avoids errors caused on
# some systems.
class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))

class CommandLine:
    """
    authors: Darrin Schultz
    Handle the command line, usage and help requests.
    """

    def __init__(self) :
        """For stage 1, the necessary arguments for the MerParse class are:
          - inputFile
          - sweep
          - sList
          - asPrefix
          - asSI
          - genus
          - species
        """

        self.parser=argparse.ArgumentParser(description=__doc__)
        self.parser.add_argument("-i", "--inputConfig",
                            type=str,
                            action=FullPaths,
                            required=True,
                            help="""The meraculous config file upon which the
                            sweep will be based.""")
        self.parser.add_argument("-g", "--genus",
                            type=str,
                            help="""The genus name for the sample that will be used
                            for naming config files and directory names.""")
        self.parser.add_argument("-s", "--species",
                            type=str,
                            help="""The species name for the sample that will be used
                            for naming config files and directory names.""")
    def parse(self):
        self.args = self.parser.parse_args()
        print(self.args)

def main():
    """
    1. Reads in a meraculous config file and outputs glotk project files
    """
    parser = CommandLine()
    #this block from here: http://stackoverflow.com/a/4042861/5843327
    if len(sys.argv)==1:
        parser.parser.print_help()
        sys.exit(1)
    parser.parse()
    myArgs = parser.args
    print(myArgs)

    # 1. Verify that the current directory isn't already a gloTK project
    cwd = os.path.abspath(os.getcwd())
    if gloTK.utils.dir_is_glotk(cwd):
        raise ValueError("""ERROR: Are you in an empty directory? This appears
        to be a gloTK project directory already. Please try again in a new
        directory!""")

    # 2. Reads in a meraculous config file and outputs all of the associated config
    #    files to $PWD/glotk_info
    configFile = ConfigParse(myArgs.inputConfig)
    gloTK_info=os.path.join(cwd, "gloTK_info")
    gloTK.utils.safe_mkdir(gloTK_info)
    shutil.copyfile(myArgs.inputConfig, os.path.join(gloTK_info, "project_init.config"))
    configFile.save_yaml(os.path.join(gloTK_info, "input_config.yaml"))

    # 3. Make symlinks for each read pair
    gloTK_reads = os.path.join(cwd, "gloTK_reads")
    gloTK.utils.safe_mkdir(gloTK_reads)
    reads0 = os.path.join(gloTK_reads, "reads0")
    gloTK.utils.safe_mkdir(reads0)
      # the files get saved in `project_dir/glotk_reads/reads0
    params_new = configFile.sym_reads_new_config(reads0, sym=True)
      # the yaml config files get saved in `project_dir/glotk_info/read_configs/`
    read_params = os.path.join(gloTK_info, "read_configs")
    gloTK.utils.safe_mkdir(read_params)
    params_new = configFile.save_yaml(os.path.join(read_params, "reads0.yaml"))

if __name__ == "__main__":
    sys.exit(main())
