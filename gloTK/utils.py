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
  - defines short functions to perform a general task

"""
import subprocess
import os
import gzip
from Bio import SeqIO
from Bio.SeqUtils import GC
from collections import Counter
import warnings

def fastq_info(path):
    """ Found some info about how to ignore warnings in code blocks here:
      - http://stackoverflow.com/questions/14463277/how-to-disable-python-warnings
    """
    numBases = 0
    numReads = 0
    readLengths = Counter()
    GCTot = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        handle = gzip.open(path, "rt")
        for record in SeqIO.parse(handle, "fastq"):
            numBases += len(record)
            numReads += 1
            readLengths[len(record)] += 1
            GCTot += sum(record.seq.count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
        handle.close()
    GCPer = (GCTot/numBases)
    avgReadLen = (sum(value*count for value,count in readLengths.items())/numReads)
    return {"numBases": numBases,
            "numReads": numReads,
            "numGCBases": GCTot,
            "portionGC": GCPer,
            "avgReadLen": avgReadLen}
