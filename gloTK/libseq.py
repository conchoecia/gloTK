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
import glob

class LibSeq(UserDict):
    def __init__(self, line):
        self.indices = {
                  1: "wildcard",
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
        self.data = {self.indices.get(x): "" for x in self.indices}
        self.libseq_parse(line)

    def libseq_parse(self, line):
        """This method parses a line containing a lib_seq string and adds it to
        the lib_seq key in the self.params dictionary"""
        count = line.count(',')
        if count > 1:
            raise ValueError("""ERROR: There are too many commas in your
            lib_seq line. There should only be one, between the two filepath
            globs.""")
        if count == 0:
            raise ValueError("""ERROR: There are no commas in your
            lib_seq line. There should be one, between the two filepath
            globs.""")
        line = [x.strip() for x in line.split() if x]

        #This block tries to fix user input if there is an inappropriate
        # number of commas or spaces associated with commas
        if len(line) > 13:
            commas = [x for x in line if "," in x]
            if commas:
                #this handles the case if there was an extra space added on
                # either side of the comma
                if (len(commas) == 1) and (commas[0] != ','):
                    if (commas[0][-1] == ',') or (commas[0][0] == ','):
                        line = [line[0]] + ["".join(line[1:3])] + line[3:]
                #This handles the case where there is an extra space on either
                #side of the comma
                elif line[2] == ',':
                    line = [line[0]] + ["".join(line[1:4])] + line[4::]

        #make sure there is an appropriate numer of elements in list
        if (len(line) ) != 13:
            print(line)
            print(commas)
            raise ValueError("""ERROR: There are {0} values in the lib_seq line when
                             there should be 12. Check your input meraculous
                             config file for errors and make sure that it matches
                             the specification for the Meraculous manual.""".format(
                                 len(line)))
        for index in self.indices:
            if index == 1:
                globs = [x.strip() for x in line[1].split(',') if x]
                #make sure there are only two filepaths to search for globs
                if len(globs) != 2:
                    raise ValueError("""ERROR: remember to split your glob
                    filepaths in lib_seq with a comma. This error occured because
                    the program was looking for two glob filepaths, but instead
                    found {}: {}""".format(len(globs), globs))
                self.data["globs"] = globs
                forwards = glob.glob(globs[0])
                reverses = glob.glob(globs[1])
                if not len(forwards) == len(reverses):
                    raise ValueError("""ERROR: the number of reverse read files does not
                    match the number of forward read files""")
                self.data["pairs"] = [x for x in zip(sorted(forwards), sorted(reverses))]
            self.data[self.indices.get(index)] = line[index]

    def __repr__(self):
        print_str = ""
        for index in sorted(self.indices):
            value = self.data.get(self.indices.get(index))
            print_str += "{0} ".format(value.replace(" ",""))
        return print_str[0:-1]
