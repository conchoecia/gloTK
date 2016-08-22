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


    def test_CommandLine_sweep_known(self):
        """Try to perform sweep on all the options.
        Namespace for sweep should contain the option.
        """
        args = ["-i", self.configPath, "--sweep", "mer_size",
                "--sstart", "1",
                "--sstop", "2",
                "--sinterval", "1"]
        parsedArgs = CommandLine(args).args
        self.assertEqual("mer_size", parsedArgs.sweep)


if __name__ == '__main__':
    unittest.main()
