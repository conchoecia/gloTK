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
from gloTK import Seqtk

import os
import subprocess
import gzip

class seqtk_test_case(unittest.TestCase):
    """Tests that the seqtk class works correctly"""
    def setUp(self):
        self.readPath = os.path.join(os.path.abspath(os.path.dirname(__file__)),"phix174Test/reads/")

    def test_sample(self):

        """This method tests that seqtk works correctly by subsampling a small
        test dataset. It creates a zipped file and then verifies that
        the unzipped file has the correct number of lines
        """
        self.inpath = os.path.join(self.readPath, "SRR353630_2500_1.fastq.gz")
        self.outpath = os.path.join(self.readPath, "sampled.txt.gz")
        sampled = Seqtk(100, self.inpath, 2, self.outpath)
        sampled.sample()
        linecount = 0
        with gzip.open(self.outpath, 'rb') as f:
             for i, l in enumerate(f):
                 pass
             linecount = i + 1
        self.assertEqual(linecount, 8)
        os.remove(self.outpath)


if __name__ == '__main__':
    unittest.main()
