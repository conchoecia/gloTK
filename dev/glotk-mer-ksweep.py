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

"""

title: glotk-mer-ksweep.py
authr: darrin schultz
vrsin: 0.01

v0.01 - Tue Aug  9 18:28:25 PDT 2016
      - Started file

This program will do a few things:
1. It will read in a meraculous config file and parse the input.
2. It will write new config files from the input and perform a multi-assembly
   parameter sweep. Simply put, it runs many assemblies through a user-defined
   parameter range so that the user can find the optimum assembly.
"""
from collections import UserDict

class libSeq(UserDict):
    def __init__(self):
        values = ["wildcard", "name", "insertAvg", "insertSdev",
                  "hasInnieArtifact", "isRevComped", "useForContiging",
                  "scaffRound", "useForGapClosing",
                  "5p_wiggleRoom", "3p_wiggleRoom"]
        self.data = {x: "" for x in values}

    def __repr__(self):
        print_str = "{"
        for key in self.data:
            val = self.data.get(key)
            print_str += "'{0}': '{1}',".format(key, val)
        print_str = print_str[:-1]
        print_str += "}"
        return print_str

class merParse:
    """This class:

    1. reads in an input file and extrapolates all of the possible
       parameters for Meraculous
    2. 

    """
    def __init__(inputFile):
        self.input = inputFile
        self.params = None
