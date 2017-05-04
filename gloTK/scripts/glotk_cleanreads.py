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

"""

title: glotk-cleanreads
authr: darrin schultz

This program:
1. Takes a read set from a gloTK project and modifies the reads using SeqPrep or
   Trimmomatic
2. Outputs a new directory with the reads in the glotk_reads/ directory
3. Outputs a yaml file for those reads in the glotk_info/read_configs/ directory
4. Logs the activity in the glotk_info/activity_log/ directory (??)
5. Generates a fastqc report for the reads

Usage:
 < glotk-project --readsNum 0 --outNum 1 SeqPrep2 >
   - runs SeqPrep2 with the defaults,o
     * -q 13
     * -L 30
     * -A AGATCGGAAGAGCACACGTC
     * -B AGATCGGAAGAGCGTCGTGT
     * -d 1
     * -C AGATCGGAAGAGCACACGTC
     * -D AGATCGGAAGAGCGTCGTGT
 < glotk-project --readsNum 0 --outNum 1 SeqPrep2 -A GATTACA -B TGTAATC >
   - runs SeqPrep2 with modified forward and reverse strands

Notes:
- Currently does not support merged output

"""

#import things for rest of program
import argparse
import copy
import os
import shutil
import sys

#multiprocessing stuff
from functools import partial
from multiprocessing import cpu_count
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

#import gloTK stuff
from gloTK import ConfigParse
import gloTK.utils
import gloTK.wrappers

#This class is used in argparse to expand the ~. This avoids errors caused on
# some systems.
class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))

class CommandLine:
    """
    authors: Darrin Schultz
    Handle the command line, usage and help requests.
    """

    def __init__(self) :
        self.parser=argparse.ArgumentParser(description=__doc__,
                                            formatter_class=argparse.RawDescriptionHelpFormatter)
        self.parser.add_argument("-r", "--readsNum",
                            type=int,
                            required=True,
                            help="""Use the reads at this location for the
                            sequence processing""")
        self.parser.add_argument("-n", "--outNum",
                            type=int,
                            required=True,
                            help="""The processed reads will be output to this
                            location""")
        subparsers = self.parser.add_subparsers(title = "Sequence Cleaner Operations",
                            description = "These are arguments for SeqPrep2",
                            #the dest option passes the subparser name to this
                            dest="operation",
                            help='sub-command help')

        SPParser = subparsers.add_parser("SeqPrep2", help='SeqPrep2 help')
        SPParser.add_argument("-q", "--qualCutoff",
                              type = int, default = 13,
                              help = """quality score cutoff for mismatches to
                              be counted in overlap; default = 13""")
        SPParser.add_argument("-L", "--lenCutoff",
                              type = int, default = 30,
                              help = """Minimum length of a trimmed or merged
                              read to print it; default = 30""")
        SPParser.add_argument("-A", "--forAdapter",
                              type = str, default = "AGATCGGAAGAGCACACGTC",
                              help = """forward read primer/adapter sequence
                              to trim as it would appear at the end of
                              a read (recommend about 20bp of this)
                              (should validate by grepping a file);
                              default (genomic non-multiplexed
                              adapter1) = AGATCGGAAGAGCACACGTC""")
        SPParser.add_argument("-B", "--revAdapter",
                              type = str, default = "AGATCGGAAGAGCGTCGTGT",
                              help = """reverse read primer/adapter sequence to
                              trim as it would appear at
                              the end of a read (recommend about 20bp
                              of this) (should validate by grepping a
                              file); default (genomic non-multiplexed
                              adapter2) = AGATCGGAAGAGCGTCGTGT""")
        SPParser.add_argument("-d", "--editReject",
                              type = int, default = 1,
                              help = """perform sequence match rejection and set
                              cutoff edit distance for rejection sequences;
                              default = 1""")
        SPParser.add_argument("-C", "--forwardReject",
                              type = str, default = "AGATCGGAAGAGCACACGTC",
                              help = """first read primer rejection sequence;
                              default = AGATCGGAAGAGCACACGTC""")
        SPParser.add_argument("-D", "--reverseReject",
                              type = str, default = "AGATCGGAAGAGCGTCGTGT",
                              help = """second read primer rejection sequence;
                              default = AGATCGGAAGAGCGTCGTGT""")
        TMParser = subparsers.add_parser("trimmomaticSE", help='trimmomatic help')
        # TMParser.add_argument("--ILLUMINACLIP",
        #                       help = """cut adapter and other illumina-specific
        #                       seqs.""")
        TMParser.add_argument("--SLIDINGWINDOW",
                              help = """perform a sliding window trimming, cutting
                              once the avg quality within the window falls below a
                              threshold.""")
        TMParser.add_argument("--LEADING",
                              help = """cut bases off the start of a read,
                              if below a threshold quality""")
        TMParser.add_argument("--TRAILING",
                              help = """cut bases off the end of a read, if below
                              a threshold quality""")
        TMParser.add_argument("--CROP",
                              help = """cut the read to a specified length""")
        TMParser.add_argument("--HEADCROP",
                              help = """cut the specified number of bases from the
                              start of the read""")
        TMParser.add_argument("--MINLEN",
                              help = """drop the read if it is below a specified
                              length""")
        TMParser.add_argument("--jarPath",
                              type = str,
                              required = True,
                              help = """input the path to the Trimmomatic jar
                              file""")
        PEParser = subparsers.add_parser("trimmomaticPE", help='trimmomatic help')
        PEParser.add_argument("--SLIDINGWINDOW",
                              help = """perform a sliding window trimming, cutting
                              once the avg quality within the window falls below a
                              threshold.""")
        PEParser.add_argument("--LEADING",
                              help = """cut bases off the start of a read,
                              if below a threshold quality""")
        PEParser.add_argument("--TRAILING",
                              help = """cut bases off the end of a read, if below
                              a threshold quality""")
        PEParser.add_argument("--CROP",
                              help = """cut the read to a specified length""")
        PEParser.add_argument("--HEADCROP",
                              help = """cut the specified number of bases from the
                              start of the read""")
        PEParser.add_argument("--MINLEN",
                              help = """drop the read if it is below a specified
                              length""")
        PEParser.add_argument("--jarPath",
                              type = str,
                              required = True,
                              help = """input the path to the Trimmomatic jar
                              file""")






    def parse(self):
        self.args = self.parser.parse_args()
        # Verify that you have selected a valid option for readcleaner
        cleanreadsOptions = ["SeqPrep2", "trimmomaticSE"]
          # The operation is the first item in the string, like SeqPrep2 or Trimmomatic
        selectedOperation = self.args.operation
          # The operation options are what to pass to the operation
        if not selectedOperation in cleanreadsOptions:
            raise ValueError("""The operation you selected, {}, is not one of the
            options, {}""".format(selectedOperation, cleanreadsOptions))

        # Verify that we are in a gloTK directory
        setattr(self.args, "gloTKDir", os.getcwd())
        if not gloTK.utils.dir_is_glotk(self.args.gloTKDir):
            raise ValueError("""This is not a gloTK directory. Please initialize
            this directory with glotk-project.""")

def main():
    # Parse the arguments
    parser = CommandLine()

    #this block from here: http://stackoverflow.com/a/4042861/5843327
    if len(sys.argv)==1:
        parser.parser.print_help()
        sys.exit(1)
    parser.parse()
    options = parser.args
    print(options)

    # Verify that the reads that you have chosen exist
    exists = gloTK.utils.reads_and_yaml_exist(options.gloTKDir, options.readsNum)
    if not exists[0]:
        raise ValueError(exists[1])

    # Load the yaml file for the reads
    yamlPath = os.path.join(options.gloTKDir,
                            "gloTK_info/read_configs/reads{}.yaml".format(
                                options.readsNum))
    readsYaml = ConfigParse(yamlPath)

    # Check that the new read directory and yaml file don't already exist
    exists = gloTK.utils.reads_and_yaml_exist(options.gloTKDir, options.outNum)
    if exists[0]:
        raise ValueError(exists[1])

    # Make a directory for the new processed reads
    newReadsDir = os.path.join(options.gloTKDir, "gloTK_reads/reads{}".format(
                                                               options.outNum))
    #Set the outdir attribute for log file
    setattr(options, "outdir", newReadsDir)

    #make a list of instances to run in parallel
    instances = []
    #Make the directory to write the new files to
    gloTK.utils.safe_mkdir(newReadsDir)
    for libseq in readsYaml.params["lib_seq"]:
        for pair in libseq["pairs"]:
            pairArgs = vars(options)
            pairArgs["forwardPath"] = pair[0]
            pairArgs["reversePath"] = pair[1]
            pairArgs["forwardOutFile"] = os.path.join(
                newReadsDir, os.path.basename(pair[0]))
            pairArgs["reverseOutFile"] = os.path.join(
                newReadsDir, os.path.basename(pair[1]))
            if options.operation == "SeqPrep2":
                #The ** turns it into kwargs http://stackoverflow.com/questions/1559638/
                thisInstance = gloTK.wrappers.Seqprep(**pairArgs)
                instances.append(thisInstance)
            if options.operation == "trimmomaticSE":
                #there are two instances to append here since the software
                # only operates on one pair at a time
                pairArgs["input"] = pairArgs["forwardPath"]
                pairArgs["output"] = pairArgs["forwardOutFile"]
                thisInstance = gloTK.wrappers.TrimmomaticSE(**pairArgs)
                instances.append(thisInstance)
                pairArgs["input"] = pairArgs["reversePath"]
                pairArgs["output"] = pairArgs["reverseOutFile"]
                thisInstance = gloTK.wrappers.TrimmomaticSE(**pairArgs)
                instances.append(thisInstance)
            if options.operation == "trimmomaticPE":
                #there are some special circumstances for PairedEnd Data
                pairArgs["forwardOutFilePaired"] = pairArgs["forwardOutFile"]
                pairArgs["forwardOutFileUnpaired"] = pairArgs["reverseOutFile"]
                pairArgs["reverseOutFilePaired"] = pairArgs["forwardOutFile"].replace(".fastq.gz", ".unpaired.fastq.gz")
                pairArgs["reverseOutFileUnpaired"] = pairArgs["reverseOutFile"].replace(".fastq.gz",".unpaired.fastq.gz")
                thisInstance = gloTK.wrappers.TrimmomaticPE(**pairArgs)
                instances.append(thisInstance)

    #process the reads
    print("Using {} processors.".format(len(instances)))
    # run the program for each instance
    # pool size is the number of simultaneous runs for the server
    threadpoolsize = len(instances)
    if cpu_count() < threadpoolsize:
        threadpoolsize = cpu_count()
    print(threadpoolsize)
    pool = ThreadPool(threadpoolsize)
    print(len(instances))
    results = pool.map(seqProcess_run_helper, instances)
    pool.close()
    pool.join()

    #Make a new ConfigParse object to modify the glob and reads
    newYaml = readsYaml.sym_reads_new_config(newReadsDir, sym=False, mv=False)
    newYaml.save_yaml(os.path.join(options.gloTKDir,
                            "gloTK_info/read_configs/reads{}.yaml".format(
                                options.outNum)))

    # run fastqc on all of the files
    #set the outdir for the fastqc files and make that dir
    fastqcOutDir = os.path.join(options.gloTKDir,
                                "gloTK_fastqc/reads{}".format(options.outNum))
    #fastqc doesn't make the directory if it doesn't exist
    gloTK.utils.safe_mkdir(fastqcOutDir)
    #turn the read list into a string
    readlist = newYaml.all_reads()
    print("readlist: ", readlist)
    #determine threads to use
    threadpoolsize = len(readlist)
    if cpu_count() < threadpoolsize:
        threadpoolsize = cpu_count()
    #run fastqc on all the reads
    args = {"readlist": readlist,
            "outdir": fastqcOutDir,
            "threads": threadpoolsize}
    gloTK.wrappers.Fastqc(**args)

def seqProcess_run_helper(instance):
    """This method is a helper method for class MerRunner. It allows
    multiprocessing module in python to map parallelism to a class method."""
    instance.run()

if __name__ == "__main__":
    sys.exit(main())
