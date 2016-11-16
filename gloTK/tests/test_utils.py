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
This class tests the classes and methods for utils.py

Thanks to Mark Howison and Casey Dunn @ Brown for inspiration and direct source
for some of the code in this file. See Biolite @ https://bitbucket.org/caseywdunn/biolite
"""

import unittest
from gloTK.wrappers import Seqprep
import gloTK.utils

import gzip
import os
import shutil
import subprocess

class utils_test_case(unittest.TestCase):
    """Tests that the functions in utils work correctly"""
    def setUp(self):
        self.readPath = os.path.join(os.path.abspath(os.path.dirname(__file__)),"phix174Test/reads/")
        self.forwardPath = os.path.join(self.readPath, "SRR353630_2500_1.fastq.gz")

    def test_fastqParse(self):
        result = gloTK.utils.fastx_basename(self.forwardPath)
        self.assertEqual("SRR353630_2500_1", result)

    def test_unmerged(self):
        """ verify that fastq_info is working correctly
           ~Correct values~
            "numBases": 375000,
            "numReads": 2500,
            "numGCBases": 171574,
            "avgReadLen": 150}
        """
        #correct values
        correct = {"numBases": 375000,
                   "numReads": 2500,
                   "numGCBases": 171574,
                   "avgReadLen": 150}

        info = gloTK.utils.fastq_info(self.forwardPath)

        for key in correct:
            self.assertEqual(correct[key], info[key])

    def test_is_gloTK(self):
        """Tests that gloTK.utils.dir_is_glotk() works by first verifying in the
        negative case, then the positive case. This method cleans up the
        directories it creates."""
        testcwd = os.path.dirname(os.path.abspath(__file__))
        os.chdir(testcwd)
        testgloTK = os.path.join(testcwd, "gloTK_assemblies")
        if os.path.isdir(testgloTK):
            shutil.rmtree(testgloTK)
        self.assertEqual(False, gloTK.utils.dir_is_glotk(testcwd))
        gloTK.utils.safe_mkdir(testgloTK)
        #IMPORTANT - do not chage testgloTK in shutil.rmtree to another dir
        self.assertEqual(True, gloTK.utils.dir_is_glotk(testgloTK))
        shutil.rmtree(testgloTK)

if __name__ == '__main__':
    unittest.main()
