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

class Seqtk:
    def __init__(self, seed, inputfile, readcount, outputfile):
        self.seed = seed
        self.inputfile = inputfile
        self.readcount = readcount
        self.outputfile = outputfile

    def sample(self):
        """Run seqtk to sample the reads and return output if there is any"""
        try:
            os.remove(self.outputfile)
        except OSError:
            pass

        call = ['seqtk', 'sample',
                '-s', "{}".format(self.seed),
                self.inputfile, "{}".format(self.readcount)]
        call_this = " ".join(call)
        sample = subprocess.run(call_this,
                                shell=True,
                                stdout=subprocess.PIPE).stdout
        with gzip.open(self.outputfile, 'wb') as f:
            f.write(sample)

        # compress = subprocess.Popen("gzip -c > {}".format(self.outputfile),
        #                           shell=True,
        #                           stdin=subprocess.PIPE, stdout=subprocess.PIPE
        #                           )
        # compress.stdin.write(compress.stdout)
        # compress.communicate()

        #.stdout.decode("utf-8").split("\n").strip()

        #return compress
