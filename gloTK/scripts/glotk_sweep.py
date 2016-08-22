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

#This class is used in argparse to expand the ~. This avoids errors caused on
# some systems.
class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))

class CommandLine:
    """
    authors: David Bernick and Darrin Schultz
    Handle the command line, usage and help requests.
    """
    def __init__(self) :
        """For stage 1, the necessary arguments for the MerParse class are:
          - inputFile
          - sweep
          - sStart
          - sStop
          - sInterval
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
        self.parser.add_argument("--sstart",
                            help="""The start for the sweep parameter""")
        self.parser.add_argument("--sstop",
                            help="""The stop for the sweep parameter""")
        self.parser.add_argument("--sinterval",
                            help="""The interval between the sweep start and stop.""")
        self.parser.add_argument("-p", "--prefix",
                            type=str,
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


        self.args = self.parser.parse_args()

        # check for input errors
        self._check_errors()

    def _check_errors(self):
        """Call this function after parsing the args to see if there are any
        errors in the way things are input. Specifically for glotk-sweep, make
        sure that all of the parameters for sweep have arguemnts if at least one
        does."""
        # use this to make sure that sweep, sstart, sstop, and sinterval
        # all have values if one of them has values
        # (P^Q^R^S) OR (-P and -Q and -R and -S)
        # above is the minimal form for this logical statement in DNF
        if not ((self._inn(self.args.sweep) and self._inn(self.args.sstart)
                 and self._inn(self.args.sstop) and self._inn(self.args.sinterval))
                or
                (not self._inn(self.args.sweep) and not self._inn(self.args.sstart) and not
                 self._inn(self.args.sstop) and not self._inn(self.args.sinterval))):
            mutually_inclusive = ["sweep", "sstart", "sstop", "sinterval"]
            print_str = ""
            for each in mutually_inclusive:
                print_str += "    {0}: {1}\n".format(each, getattr(self.args, each))
            raise AssertionError("""You specified one or more of the --sweep,
            --sstart, --stop, or --sinterval arguments but did not specify all of
            them. All of them are required when running the program in sweep mode.
            Do not specify any of these arguments if not running the program in
            sweep mode.\n{0}""".format(print_str))

    def _inn(self, value):
        """Checks to verify that an argument is not None, aka that it has a value.
        Lets you use boolean logic on arguments since Python doesn't have
        mutually inclusive groups."""
        return value is not None

def mer_runner_dummy(instance):
    """This method is a helper method for class MerRunner. It allows
    multiprocessing module in python to map parallelism to a class method."""
    instance.meraculous_runner()

class MerRunner:
    """This class has one instance per Meraculous run and is accessed with the
    partial module"""
    def __init__(self, runName, configPath, cleanup=0):
        """The cleanup parameter is what is passed to the run_meraculous script"""
        self.runName = runName
        self.configPath = configPath
        self.cleanup =  cleanup
        self.cwd = os.path.abspath(os.getcwd())
        self.allAssembliesDir = os.path.join(self.cwd, "assemblies")
        self.thisAssemblyDir = os.path.join(self.allAssembliesDir, self.runName)
        self.reportsDir = os.path.join(self.cwd, "reports")

        self.callString = "bash run_meraculous.sh -c {0} -dir {1} -cleanup_level {2}".format(
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
    print(sys.argv)
    if len(sys.argv) <= 1:
        print("""Please input some options for this program or see how it is used.""")
        parser = CommandLine(["--help"])
    else :
        parser = CommandLine()
        myArgs = parser.args
        print(myArgs)

    #Figure out how many processors to give to each assembly since we will be
    # running some things in parallel. The MerParse class will handle overriding
    # whatever is found in the config file in the read_config() method.
    procsPerAssembly = int(myArgs.maxProcs / myArgs.simultaneous)
    setattr(myArgs, "maxProcs", procsPerAssembly)

    # 1. Reads in a meraculous config file and outputs all of the associated config
    #    files to $PWD/configs

    merparser = MerParse(myArgs.inputConfig, myArgs.sweep, myArgs.sstart,
                         myArgs.sstop, myArgs.sinterval, myArgs.maxProcs,
                         asPrefix = myArgs.prefix,
                         asSI = myArgs.index,
                         genus = myArgs.genus,
                         species = myArgs.species)
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
        thisInstance = MerRunner(runName, configPath)
        instances.append(thisInstance)

    if len(instances) == 0:
        print("There are no meraculous folders in this directory. Exiting")
    elif len(instances) > 0:
        # run the program for each instance
        # pool size is the number of simultaneous runs for the server
        pool = ThreadPool(myArgs.simultaneous)
        results = pool.map(mer_runner_dummy, instances)
        pool.close()
        pool.join()

if __name__ == "__main__":
    sys.exit(main())
