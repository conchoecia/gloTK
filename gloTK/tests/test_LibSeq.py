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
import ast


class libSeq_test_case(unittest.TestCase):
    """Tests that libSeq works correctly, like a dictionary
    """
    def test_libSeq_construction(self):
        """This verifies that the LibSeq class is constructed with all of the
        parameters that are present in the meraculous config files"""
        myLib = LibSeq()
        values = ["wildcard", "name", "insertAvg", "insertSdev", "avgReadLn",
                  "hasInnieArtifact", "isRevComped", "useForContiging",
                  "scaffRound", "useForGapClosing",
                  "5p_wiggleRoom", "3p_wiggleRoom"]
        data = {x: "" for x in values}
        #go through all the keys in the LibSeq object to make sure they match
        # all the params in the Meraculous manual
        for key in myLib.keys():
            data.pop(key)
        self.assertEqual(data, {})

    def test_libSeq_repr_method(self):
        """Verify that the __repr__() methods returns the correctly-formatted
        string when using print(LibSeq())"""
        myLib = LibSeq()
        values = {"wildcard": "/my/path/*.gz",
                  "name": "CAT",
                  "insertAvg": "600",
                  "insertSdev": "1",
                  "avgReadLn": "150",
                  "hasInnieArtifact": "0",
                  "isRevComped": "0",
                  "useForContiging": "1",
                  "scaffRound": "1",
                  "useForGapClosing": "1",
                  "5p_wiggleRoom": "0",
                  "3p_wiggleRoom": "0"}
        for key in values:
            myLib[key] = values[key]
        out_str = "/my/path/*.gz CAT 600 1 150 0 0 1 1 1 0 0"
        self.assertEqual(str(myLib), out_str)

    # def test_libSeq_remove_spaces(self):
    #     """Test that values have spaces stripped upon entry in self.data."""
    #     myLib = LibSeq()
    #     values = {"wildcard": "/my/path/*.gz",
    #               "name": "CAT",
    #               "insertAvg": "600",
    #               "insertSdev": "1",
    #               "avgReadLn": "150",
    #               "hasInnieArtifact": "0",
    #               "isRevComped": "0",
    #               "useForContiging": "1",
    #               "scaffRound": "1",
    #               "useForGapClosing": "1",
    #               "5p_wiggleRoom": "0",
    #               "3p_wiggleRoom": "0"}
    #     values_spaces = {"wildcard": " /my/path/*.gz ",
    #               "name": " CAT ",
    #               "insertAvg": " 600 ",
    #               "insertSdev": " 1 ",
    #               "avgReadLn": " 150 ",
    #               "hasInnieArtifact": " 0 ",
    #               "isRevComped": " 0 ",
    #               "useForContiging": " 1 ",
    #               "scaffRound": " 1 ",
    #               "useForGapClosing": " 1 ",
    #               "5p_wiggleRoom": " 0 ",
    #               "3p_wiggleRoom": " 0 "}
    #     for key in values_spaces:
    #         myLib[key] = values_spaces[key]
    #     #the dict obkect in myLib should equal the values dict above.
    #     # aka all extraneous spaces should be removed
    #     # this is a redundant feature in the software since the MerParse
    #     # object splits on spaces, but this is useful if I need to modify
    #     # LibSeq objects with other classes.
    #     self.assertEqual(myLib, values)


if __name__ == '__main__':
    unittest.main()
