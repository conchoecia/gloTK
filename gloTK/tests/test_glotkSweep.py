#!/usr/bin/env python3
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

"""@author Darrin Schultz
This class tests the classes and methods for glotk_sweep.py
"""

import unittest
import argparse
from gloTK.scripts.glotk_sweep import CommandLine
import os
import shutil

#call the assembly with the shell
import subprocess


class ErrorRaisingArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        #print(message)
        raise ValueError(message)  # reraise an error

class sweep_test_case(unittest.TestCase):
    """Tests that the merParse class works correctly"""

    def setUp(self):
        self.parser=ErrorRaisingArgumentParser()
        self.configPath = os.path.join(os.path.abspath(os.path.dirname(__file__)),"ksweep_test.config")
        self.parser.add_argument(
            "-s", "--sweep",
            type=str,
            choices=["mer_size", "bubble_depth_threshold"],
            required=True)

    def test_required_unknown(self):
        """Try to perform sweep on something that isn't an option.
        Argparse will raise an error and the test will catch it passing.
        This is more so for demonstration purposes, since the CommandLine class
        inherits the normal argparse module without the inherited error method.
        """
        args = ["--sweep", "NADA"]
        with self.assertRaises(ValueError) as cm:
            self.parser.parse_args(args)
        print('msg:',cm.exception)
        self.assertIn('invalid choice', str(cm.exception))

class assembly_test_case(unittest.TestCase):
    """Tests a swept assembly"""
    def setUp(self):
        self.testRunDir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "phix174Test")
        self.configDir = os.path.join(self.testRunDir, "configs")
        self.assemDir = os.path.join(os.path.abspath(self.testRunDir) , "assemblies")
        self.outputParentDir = os.path.join(os.path.abspath(os.path.dirname(__file__)))
        self.read1 = os.path.join(self.testRunDir, "reads/SRR353630_2500_1.fastq.gz")
        self.read2 = os.path.join(self.testRunDir, "reads/SRR353630_2500_2.fastq.gz")
        self.readLink1 = os.path.join(self.testRunDir, "assemblies/SRR353630_2500_1.fastq.gz")
        self.readLink2 = os.path.join(self.testRunDir, "assemblies/SRR353630_2500_2.fastq.gz")


    def tearDown(self):
        """This directory is called after every method is run even if the method
        produces an error.

        Here, it is used to remove any report files that are generated while
        running these tests.
        """


    def test_sweeps(self):
        """Make sure that the program runs without error"""
        os.chdir(self.testRunDir)
        #make sure symlink for reads exists
        if not os.path.exists(self.assemDir):
            os.makedirs(self.assemDir)
        os.symlink(self.read1, self.readLink1)
        os.symlink(self.read2, self.readLink2)

        callString = ["glotk-sweep",
                      "-i", "phix174.config",
                      "-s", "mer_size",
                      "--sstart", "21",
                      "--sstop", "25",
                      "--sinterval", "2",
                      "-n", "1"]
        p = subprocess.run(callString, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True)
        output = str(p.stdout)
        err = str(p.stderr)
        print(output)
        print(err)

        self.assertTrue(os.path.isfile(self.readLink1))
        self.assertTrue(os.path.isfile(self.readLink2))
        os.remove(self.readLink1)
        os.remove(self.readLink2)

        #get the names of the config files to look for those dir names
        assemblies = os.listdir(self.configDir)
        for each in assemblies:
            assemblydir=os.path.join(self.assemDir, each)
            configfile=os.path.join(self.configDir, each)
            self.assertTrue(os.path.isdir(assemblydir))
            self.assertTrue(os.path.isfile(configfile))
        #remove the assembly dir after you're sure they're all there
        shutil.rmtree(self.assemDir)
        shutil.rmtree(self.configDir)

if __name__ == '__main__':
    unittest.main()
