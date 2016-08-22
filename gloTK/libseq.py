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

title: libseq.py
authr: Darrin Schultz

This is the LibSeq class that is used to store data for the "lib_seq" information
in the Meraculous config files.
"""

from collections import UserDict

class LibSeq(UserDict):
    def __init__(self):
        self.values = {1: "wildcard",
                  2: "name",
                  3: "insertAvg",
                  4: "insertSdev",
                  5: "avgReadLn",
                  6: "hasInnieArtifact",
                  7: "isRevComped",
                  8: "useForContiging",
                  9: "scaffRound",
                  10: "useForGapClosing",
                  11: "5p_wiggleRoom",
                  12: "3p_wiggleRoom"}
        self.data = {self.values.get(x): "" for x in self.values}

    def __repr__(self):
        print_str = ""
        for index in sorted(self.values):
            value = self.data.get(self.values.get(index))
            print_str += "{0} ".format(value.replace(" ",""))
        return print_str[0:-1]
