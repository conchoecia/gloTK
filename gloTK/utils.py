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

"""title: utils.py
authr: darrin schultz

This module:
  - defines short functions to perform a general task

"""
import inspect
import gzip
import os
import subprocess
import sys
import time
import warnings

from Bio import SeqIO
from Bio.SeqUtils import GC
from collections import Counter
from traceback import print_stack


def gzip_verify(filepath):
    if filepath.strip()[-3::] != ".gz":
        return "{}.gz".format(filepath.strip())
    else:
        return filepath

def info(*messages):
    """
    Prints the current GloTK module and a `message`.
    Taken from biolite
    """
    sys.stderr.write("%s.%s: " % get_caller_info())
    sys.stderr.write(' '.join(map(str, messages)))
    sys.stderr.write('\n')

def safe_mkdir(path):
    """
    Creates the directory, including any missing parent directories, at the
    specified `path`.

    Aborts if the path points to an existing regular file.

    Returns the absolute path of the directory.
    """
    if os.path.isfile(path):
        die("'{0}' is a regular file: can't overwrite" % path)
    elif os.path.isdir(path):
        info("directory '%s' already exists" % path)
    else:
        info("creating directory '%s'" % path)
        try:
            os.makedirs(path)
        except OSError as e:
            die("""failed to recursively create the directory
%s
%s
  Do you have write permision to that path?
  Or does part of that path already exist as a regular file?""" % (path, e))
    return os.path.abspath(path)

def dir_is_glotk(path):
    """check that the current directory is a glotk project folder"""
    test_set = set(["gloTK_info", "gloTK_assemblies",
                    "gloTK_configs", "gloTK_reads",
                    "gloTK_fastqc", "gloTK_kmer",
                    "gloTK_reports"])
    #http://stackoverflow.com/questions/11968976/
    files = set([f for f in os.listdir('.') if os.path.isdir(f)])
    intersection = test_set.intersection(files)
    if len(intersection) > 0:
        return True
    else:
        return False

def timestamp():
    """
    Returns the current time in :samp:`YYYY-MM-DD HH:MM:SS` format.
    """
    return time.strftime("%Y-%m-%d %H:%M:%S")

def get_caller_info(depth=2, trace=False):
    """
    Uses the inspect module to determine the name of the calling function and
    its module.

    Returns a 2-tuple with the module name and the function name.

    From Biolite package
    """
    try:
        frame = inspect.stack()[depth]
    except:
        die("could not access the caller's frame at stack index %d" % depth)
    if trace:
        print_stack(frame[0].f_back)
    func = frame[3]
    module = inspect.getmodule(frame[0])
    if module:
        return (module.__name__, func)
    else:
        return ('<unknown>', func)

def die(*messages):
    """
    Prints the current BioLite module and an error `message`, then aborts.
    """
    sys.stderr.write("%s.%s: " % get_caller_info(trace=True))
    sys.stderr.write(' '.join(map(str, messages)))
    sys.stderr.write('\n')
    sys.exit(1)


def fastx_basename(path):
    split = os.path.splitext(os.path.basename(path))
    noZone = [".fastq",".fq",".fasta", ".fa",
               ".gz",".gzip", ".gzipped"]
    while split[1] in noZone:
        split = os.path.splitext(os.path.basename(split[0]))

    return split[0]

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
