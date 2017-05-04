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
from gloTK import ConfigParse


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
        #This printed out too many times for anyone's good.
        #info("directory '%s' already exists" % path)
        pass
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

def fastq_append(filename, appendThis):
    """This removes the ".fastq.gz", or other ".A.B" from a file,
    appends something "_appendThis" to the file name, then replaces the ".A.B"
    """
    dirname = os.path.dirname(filename)
    basename = os.path.basename(filename)
    split = basename.split(".")
    #There are two circumstances here. One is where this is a fastq.gz and one
    # where the file is just a fastq. Determine which then proceed
    if split[-1] == "gz":
        gzipped = True
        fileType = split[-2]
    else:
        gzipped = False
        fileType = split[-1]
    if fileType not in ["fastq", "fq"]:
        raise ValueError("This is not a fastq file: {}".format(basename))
    if gzipped:
        return os.path.join(dirname, "{}_{}.{}".format("".join(split[:-2]), appendThis, ".".join(split[-2:])))
    else:
        return os.path.join(dirname, "{}_{}.{}".format("".join(split[:-1]), appendThis, split[-1]))


def reads_and_yaml_exist(gloTKDir, readNum):
    """check that a yaml file and its corresponding reads exist.
    Returns a tuple of (boolean, 'return message') """
    yamlFile = os.path.join(gloTKDir, "gloTK_info/read_configs/reads{}.yaml".format(readNum))
    #if the yaml file doesn't exist, then it is false that reads and yaml exist
    if not os.path.exists(yamlFile):
        return (False, "The yaml file for reads{} does not exist.".format(yamlFile))
    thisYaml = ConfigParse(yamlFile)
    for libseq in thisYaml.params["lib_seq"]:
        for pair in libseq["pairs"]:
            for read in pair:
                if not os.path.exists(read):
                    return(False, "The read file {} does not exist".format(read))
    return (True, "Both yaml file and corresponding reads exist.")

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
    """This needs a description"""
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
