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

title: glotk-sweep
authr: darrin schultz

This program:
1. Reads in a meraculous config file and outputs all of the associated config
   files to $PWD/configs
2. The name of each run and the path to the directory is passed to a
   multiprocessing core that controls which assemblies are executed and when.
3. Each assembly is executed.

Usage: 
"--slist 21 23 57 73" to perform assemblies for kmer sizes 21, 23, 57, 73, et cetera
"""

#import things for rest of program
import sys
import argparse
import os

#stuff for call shell commands
import subprocess

#multiprocessing stuff
from functools import partial
from multiprocessing import cpu_count
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

#import gloTK stuff
from gloTK import MerParse
from gloTK import MerRunAnalyzer

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
        self.parser.add_argument("-s", "--sweep",
                            type=str,
                            choices=["mer_size", "bubble_depth_threshold"],
                            required=True,
                            help="""The parameter to sweep through in the
                            assembly.""")
        self.parser.add_argument("--slist",
                            type=str,
                            nargs='+',
                            help="""The values to sweep through""")
        self.parser.add_argument("-p", "--prefix",
                            type=str,
                            default='as',
                            help="""The assembly prefix used to name config files
                            and the output directories.""")
        self.parser.add_argument("-I", "--index",
                            type=str,
                            help="""The starting assembly index for naming the config
                            files and output directory.""")
        self.parser.add_argument("-G", "--genus",
                            type=str,
                            help="""The genus name for the sample that will be used
                            for naming config files and directory names.""")
        self.parser.add_argument("-S", "--species",
                            type=str,
                            help="""The species name for the sample that will be used
                            for naming config files and directory names.""")
        self.parser.add_argument("-n", "--simultaneous",
                            type=int,
                            default=1,
                            help="""The number of simultaneous assemblies to run.""")
        self.parser.add_argument("-M", "--maxProcs",
                            type=int,
                            default=cpu_count() - 2,
                            help="""The total number of processers used by all
                            of the assemblies combined.""")
        self.parser.add_argument("-c", "--censor",
                            type=str,
                            nargs='+',
                            help="""list of strings to censor for sending docs to
                            non-collaborators""")
        self.parser.add_argument("-q", "--quiet",
                            action='store_true')
        self.parser.add_argument("-t", "--triplet",
                            action='store_true',
                            help="""This performs diploid modes 0, 1, 2 for
                            every kmer size input""")
        self.parser.add_argument("-C", "--cleanup",
                            type=int,
                            default=1,
                            choices=[0,1,2],
                            help="""The cleanup level to pass along to the
                            run_meraculous.sh program.""")

    def parse(self):
        self.args = self.parser.parse_args()
        print(self.args)

def mer_runner_dummy(instance):
    """This method is a helper method for class MerRunner. It allows
    multiprocessing module in python to map parallelism to a class method."""
    instance.meraculous_runner()

class MerRunner:
    """This class has one instance per Meraculous run and is accessed with the
    partial module"""
    def __init__(self, runName, configPath, cleanup):
        """The cleanup parameter is what is passed to the run_meraculous script"""
        self.runName = runName
        self.configPath = configPath
        self.cleanup =  cleanup
        self.cwd = os.path.abspath(os.getcwd())
        self.allAssembliesDir = os.path.join(self.cwd, "assemblies")
        self.thisAssemblyDir = os.path.join(self.allAssembliesDir, self.runName)
        self.reportsDir = os.path.join(self.cwd, "reports")

        self.callString = "run_meraculous.sh -c {0} -dir {1} -cleanup_level {2}".format(
            self.configPath, self.runName, self.cleanup)

    def meraculous_runner(self):
        """
        Check to make sure that the allAssembliesDir has been created, if not,
        make it. This will only execute for the first time an assembly has been
        run in this directory.

        Run the directory from allAssembliesDir. The self.callString instance
        attribute tells Meraculous to name the assembly directory self.runName.

        After the run is complete, create the meraculous report, passing the
        directory containing the run (aka self.thisAssemblyDir).
        """
        #set the dir to temp assembly dir
        os.chdir(self.allAssembliesDir)

        print(self.callString)
        p = subprocess.run(self.callString, shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True)
        output = str(p.stdout)
        err = str(p.stderr)

        #generate the report for the run
        self._generate_report()

        #exit, returning the output and err
        return (output, err)

    def _generate_report(self):
        reporter = MerRunAnalyzer(self.thisAssemblyDir, self.cwd, [])
        reporter.generate_report()

def main():
    """
    1. Reads in a meraculous config file and outputs all of the associated config
       files to $PWD/configs
    2. The name of each run and the path to the directory is passed to a
       multiprocessing core that controls which assemblies are executed and when.

    """
    parser = CommandLine()
    #this block from here: http://stackoverflow.com/a/4042861/5843327
    if len(sys.argv)==1:
        parser.parser.print_help()
        sys.exit(1)
    parser.parse()
    myArgs = parser.args

    #Figure out how many processors to give to each assembly since we will be
    # running some things in parallel. The MerParse class will handle overriding
    # whatever is found in the config file in the read_config() method.
    procsPerAssembly = min(50, int(myArgs.maxProcs / myArgs.simultaneous))
    setattr(myArgs, "maxProcs", procsPerAssembly)

    # 1. Reads in a meraculous config file and outputs all of the associated config
    #    files to $PWD/configs

    merparser = MerParse(myArgs.inputConfig,
                         myArgs.sweep,
                         myArgs.slist,
                         myArgs.maxProcs,
                         asPrefix = myArgs.prefix,
                         asSI = myArgs.index,
                         genus = myArgs.genus,
                         species = myArgs.species,
                         triplet = myArgs.triplet)
    configPaths = merparser.sweeper_output()

    #make the assemblies dir ONCE to avoid a race condition for os.makedirs()
    cwd = os.path.abspath(os.getcwd())
    allAssembliesDir = os.path.join(cwd, "assemblies")
    if not os.path.exists(allAssembliesDir):
        os.makedirs(allAssembliesDir)

    #instantiate all of the classes that we will be using in parallel processing.
    # configPaths above returns a dict with the run name and abs path of config
    # as key:value pairs
    instances = []
    for runName in configPaths:
        configPath = configPaths.get(runName)
        #strip off the .config off the end of the runName, derived from configPath
        thisInstance = MerRunner(runName.strip(".config"), configPath, myArgs.cleanup)
        instances.append(thisInstance)

    if len(instances) == 0:
        print("There are no meraculous folders in this directory. Exiting")
    elif len(instances) > 0:
        print("Using {} processors.".format(procsPerAssembly * myArgs.simultaneous))
        # run the program for each instance
        # pool size is the number of simultaneous runs for the server
        pool = ThreadPool(myArgs.simultaneous)
        results = pool.map(mer_runner_dummy, instances)
        pool.close()
        pool.join()

if __name__ == "__main__":
    sys.exit(main())
