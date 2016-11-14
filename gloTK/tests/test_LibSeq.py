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
This class tests the classes and methods for glotk-mer-ksweep.py
"""

import unittest
from gloTK import LibSeq
from time import strftime as tfmt
import re
import os
import ast


class libSeq_test_case(unittest.TestCase):
    """Tests that libSeq works correctly, like a dictionary
    """
    def setUp(self):
        self.pwd = os.path.dirname(os.path.abspath(__file__))
        self.testConfig = os.path.join(self.pwd, "phix174Test/libSeq.config")
        lines=[]
        with open(self.testConfig, "r") as f:
            for line in f:
                if line.strip():
                    lines.append(line.strip())
        self.values = lines[0].split()

    def test_libSeq_construction(self):
        """This verifies that the LibSeq class is constructed with all of the
        parameters that are present in the meraculous config files"""
        dummy_line = " ".join(self.values)
        myLib = LibSeq(dummy_line)
        values = [ "wildcard", "name", "insertAvg", "insertSdev", "avgReadLn",
                  "hasInnieArtifact", "isRevComped", "useForContiging",
                  "scaffRound", "useForGapClosing",
                  "5p_wiggleRoom", "3p_wiggleRoom", "globs", "pairs"]
        data = {x: "" for x in values}
        #go through all the keys in the LibSeq object to make sure they match
        # all the params in the Meraculous manual
        for key in myLib.data.keys():
            data.pop(key)
        self.assertEqual(data, {})

    def test_wrong_num_commas(self):
        """Verifies that the wildcard string has two globs."""

        with self.assertRaises(ValueError):
            myLib = LibSeq("lib_seq phix174Test/reads/SRR353630_2500_1*.fastq.gz CAT 600 1 150 0 0 1 1 1 0 0")

    def test_libSeq_repr_method(self):
        """Verify that the __repr__() methods returns the correctly-formatted
        string when using print(LibSeq()). Also verifies that the wildcard string
        has two globs."""

        myLib = LibSeq("lib_seq phix174Test/reads/SRR353630_2500_1*.fastq.gz,phix174Test/reads/SRR353630_2500_2*.fastq.gz CAT 600 1 150 0 0 1 1 1 0 0")
        readBase=os.path.join(self.pwd, "phix174Test/reads")
        # forward=os.path.join(readBase, 
        out_str = "phix174Test/reads/SRR353630_2500_1*.fastq.gz,phix174Test/reads/SRR353630_2500_2*.fastq.gz CAT 600 1 150 0 0 1 1 1 0 0"
        print(str(myLib))
        self.assertEqual(str(myLib), out_str)

    def test_libSeq_extraspaces_right(self):
        """Verify that extra spaces won't mess up the parser"""

        myLib = LibSeq("""lib_seq  phix174Test/reads/SRR353630_2500_1*.fastq.gz ,phix174Test/reads/SRR353630_2500_2*.fastq.gz  CAT     600 1 150 0
        0 1 1 1 0 0""")
        out_str = "phix174Test/reads/SRR353630_2500_1*.fastq.gz,phix174Test/reads/SRR353630_2500_2*.fastq.gz CAT 600 1 150 0 0 1 1 1 0 0"
        self.assertEqual(str(myLib), out_str)

    def test_libSeq_extraspaces_middle(self):
        """Verify that extra spaces won't mess up the parser"""

        myLib = LibSeq("""lib_seq phix174Test/reads/SRR353630_2500_1*.fastq.gz , phix174Test/reads/SRR353630_2500_2*.fastq.gz CAT     600 1 150 0
        0 1 1 1 0 0""")
        out_str = "phix174Test/reads/SRR353630_2500_1*.fastq.gz,phix174Test/reads/SRR353630_2500_2*.fastq.gz CAT 600 1 150 0 0 1 1 1 0 0"
        self.assertEqual(str(myLib), out_str)

    def test_libSeq_extraspaces_left(self):
        """Verify that extra spaces won't mess up the parser"""

        myLib = LibSeq("""lib_seq phix174Test/reads/SRR353630_2500_1*.fastq.gz, phix174Test/reads/SRR353630_2500_2*.fastq.gz CAT     600 1 150 0
        0 1 1 1 0 0""")
        out_str = "phix174Test/reads/SRR353630_2500_1*.fastq.gz,phix174Test/reads/SRR353630_2500_2*.fastq.gz CAT 600 1 150 0 0 1 1 1 0 0"
        self.assertEqual(str(myLib), out_str)

    def test_libSeq_toomanycommas(self):
        """Verify that an error is thrown if there are too many commas"""

        with self.assertRaises(ValueError):
            myLib = LibSeq("lib_seq /my,/p,ath/*.gz,path2*.gz CAT 600 1 150 0 0 1 1 1 0 0")

    def test_libSeq_nocommas(self):
        """Verify that an error is thrown if there are no  commas"""

        with self.assertRaises(ValueError):
            myLib = LibSeq("lib_seq /my/path/*.gzpath2*.gz CAT 600 1 150 0 0 1 1 1 0 0")







if __name__ == '__main__':
    unittest.main()
