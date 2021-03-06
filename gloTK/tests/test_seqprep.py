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
This class tests the classes and methods for Seqtk in commandline_tools.py
"""

import unittest
from gloTK.wrappers import Seqprep

import os
import subprocess
import gzip

class seqprep_test_case(unittest.TestCase):
    """Tests that the seqtk class works correctly"""
    def setUp(self):
        self.readPath = os.path.join(os.path.abspath(os.path.dirname(__file__)),"phix174Test/reads/")
        self.forwardPath = os.path.join(self.readPath, "SRR353630_2500_1.fastq.gz")
        self.reversePath = os.path.join(self.readPath, "SRR353630_2500_2.fastq.gz")
        self.forwardOutFile  = ("out_forward.fastq.gz")
        self.reverseOutFile  = ("out_reverse.fastq.gz")
        self.mergedOutFile   = ("out_merged.fastq.gz")
        self.prettyOutFile   = ("pretty.txt.gz")
        self.nameVector = [self.forwardPath,
                           self.reversePath,
                           os.path.join(self.readPath, self.forwardOutFile),
                           os.path.join(self.readPath, self.reverseOutFile),
                           os.path.join(self.readPath, self.mergedOutFile),
                           os.path.join(self.readPath, self.prettyOutFile),
                           os.path.join(self.readPath, "seqprep2.log")]

    def test_unmerged(self):
        """This method tests that Seqprep is correctly processing both merged and unmerged data.
        """

        prep = Seqprep(forwardPath    = self.forwardPath,
                       reversePath    = self.reversePath,
                       forwardOutFile = self.forwardOutFile,
                       reverseOutFile = self.reverseOutFile,
                       mergedOutFile  = self.mergedOutFile,
                       prettyOutFile  = self.prettyOutFile,
                       outDir         = self.readPath)
        for path in self.nameVector:
            self.assertTrue(os.path.exists(path))
            if "SRR" not in path:
                os.remove(path)




if __name__ == '__main__':
    unittest.main()
