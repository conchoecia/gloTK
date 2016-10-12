#!/usr/bin/env python
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

"""title: commandline_tools
authr: darrin schultz

This module:
  - defines classes that interact with command line tools

Current classes
"""
import subprocess
import os
import gzip

#started 4:15PM
#ended at 5:30PM
#started at 6:25PM
# ended at 7:05PM
# started at 9PM
# ended at 12PM

class Seqprep:
    def __init__(self, forwardPath, reversePath, forwardOut, reverseOut,
                 mergedOut, prettyOut, qualCutoff = 13, lenCutoff = 30,
                 forwardAdapter ="AGATCGGAAGAGCACACGTC",
                 reverseAdapter ="AGATCGGAAGAGCGTCGTGT",
                 overlapMin = 30, editReject = 1,
                 forwardReject = "AGATCGGAAGAGCACACGTC",
                 reverseReject = "AGATCGGAAGAGCGTCGTGT"):
        self.forward       = forwardPath
        self.reverse       = reversePath
        self.forOut        = forwardOut
        self.revOut        = reverseOut
        self.mergedOut     = mergedOut
        self.prettyOut     = prettyOut
        self.prettyNum     = 50
        self.qualCutoff    = qualCutoff
        self.lenCutoff     = lenCutoff
        self.forAdapter    = forwardAdapter
        self.revAdapter    = reverseAdapter
        self.overlapMin    = overlapMin
        self.editReject    = editReject
        self.forwardReject = forwardReject
        self.reverseReject = reverseReject

    def run(self):
        call = ["SeqPrep2",
                "-f", self.forward,
                "-r", self.reverse,
                "-1", self.forOut,
                "-2", self.revOut,
                "-q", self.qualCutoff,
                "-L", self.lenCutoff,
                "-A", self.forAdapter,
                "-B", self.revAdapter,
                "-d", self.editReject,
                "-C", self.forwardReject,
                "-D", self.reverseReject]
        if self.mergedOut:
            call += ["-s", self.mergedOut,
                "-E", self.prettyOut,
                "-x", self.prettyNum,
                "-o", self.overlapMin]

        call_this = " ".join(str(x) for x in call)
        sample = subprocess.run(call_this,
                                shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)


class Seqtk:
    def __init__(self, seed, inputfile, readcount, outputfile):
        self.seed = seed
        self.inputfile = inputfile
        self.readcount = readcount
        self.outputfile = outputfile

    def run(self):
        """Run seqtk to sample the reads and return output if there is any"""
        try:
            os.remove(self.outputfile)
        except OSError:
            pass

        call = ['seqtk', 'sample',
                '-s', self.seed,
                self.inputfile, self.readcount]
        call_this = " ".join(str(x) for x in call)
        sample = subprocess.run(call_this,
                                shell=True,
                                stdout=subprocess.PIPE).stdout
        with gzip.open(self.outputfile, 'wb') as f:
            f.write(sample)
