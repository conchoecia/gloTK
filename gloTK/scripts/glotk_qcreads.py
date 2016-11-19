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

title: glotk-cleanreads
authr: darrin schultz

This program:
1. Performs one or several qc operations on a set of files.
"""

#import things for rest of program
import argparse
import copy
import os
import shutil
import sys

#multiprocessing stuff
from functools import partial
from multiprocessing import cpu_count
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

#import gloTK stuff
from gloTK import ConfigParse
import gloTK.utils
import gloTK.wrappers

class CommandLine:
    """
    authors: Darrin Schultz
    Handle the command line, usage and help requests.
    """

    def __init__(self) :
        self.parser=argparse.ArgumentParser(description=__doc__,
                                            formatter_class=argparse.RawDescriptionHelpFormatter)
        self.parser.add_argument("-r", "--readsNum",
                            type=int,
                            required=True,
                            help="""Use the reads at this location for the
                            sequence processing""")

        self.parser.add_argument("-o", "--operations",
                            nargs = '+',
                            required=True,
                            help="""This performs QC on several types of reads""")


    def parse(self):
        self.args = self.parser.parse_args()
        # Verify that you have selected a valid option for readcleaner
        qcOptions = ["fastqc", "kmergenie"]

        for operation in self.args.operations:
            if operation not in qcOptions:
                raise ValueError("""The operation you selected, {}, is not one of the
                options, {}""".format(operation, qcOptions))

        # Verify that we are in a gloTK directory
        setattr(self.args, "gloTKDir", os.getcwd())
        if not gloTK.utils.dir_is_glotk(self.args.gloTKDir):
            raise ValueError("""This is not a gloTK directory. Please initialize
            this directory with glotk-project.""")

def run_fastqc():
    readlist = inputYaml.all_reads()
    #determine threads to use
    threadpoolsize = len(readlist)
    if cpu_count() < threadpoolsize:
        threadpoolsize = cpu_count()
    #run fastqc on all the reads
    args = {"readlist": readlist,
            "outdir": operationOutDir,
            "threads": threadpoolsize}
    gloTK.wrappers.Fastqc(**args)

def run_kmergenie():
    readlist = inputYaml.all_reads()

    #make a temp directory to output the tempreads to
    tempDir = os.path.join(operationOutDir, "temp")
    gloTK.utils.safe_mkdir(tempDir)
    os.chdir(tempDir)
    tempReadFile = os.path.join(tempDir, "tempReads.fastq.gz")
    #make a concatenated read file
    with open(tempReadFile, 'wb') as outfile:
        for fname in readlist:
            with open(fname, 'rb') as infile:
                for line in infile:
                    outfile.write(line)
    os.chdir(operationOutDir)
    threadpoolsize = cpu_count() - 2
    if threadpoolsize < 1:
        threadpoolsize = 1
    #run kmergenie on the joined reads
    args = {"readFile": tempReadFile,
            "outdir": operationOutDir,
            "prefix": "reads{}".format(options.readsNum),
            "threads": threadpoolsize}
    gloTK.wrappers.Kmergenie(**args)
    # now remove the temporary directory
    os.chdir(options.gloTKDir)
    shutil.rmtree(tempDir)

def main():
    # Parse the arguments
    parser = CommandLine()

    #this block from here: http://stackoverflow.com/a/4042861/5843327
    if len(sys.argv)==1:
        parser.parser.print_help()
        sys.exit(1)
    parser.parse()
    global options
    options = parser.args
    print(options)

    yamlFilePath = os.path.join(options.gloTKDir,
                                "gloTK_info/read_configs/reads{}.yaml".format(
                                                               options.readsNum))
    global inputYaml
    inputYaml = ConfigParse(yamlFilePath)

    funcDict = {"fastqc": run_fastqc, "kmergenie": run_kmergenie}
    for operation in options.operations:
        #set the outdir for the fastqc files and make that dir
        global operationOutDir
        operationOutDir = os.path.join(options.gloTKDir,
                    "gloTK_{}/reads{}".format(operation, options.readsNum))
        # make the directory if it doesn't exist
        gloTK.utils.safe_mkdir(operationOutDir)
        #run the operation
        funcDict[operation]()

if __name__ == "__main__":
    sys.exit(main())
