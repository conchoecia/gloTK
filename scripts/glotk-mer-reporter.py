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

title: glotk-mer-reporter.py
authr: darrin schultz
vrsin: 0.10

v0.10 - Tue Aug  9 17:11:26 PDT 2016
      - Outputs results in cwd using py-gfm package
v0.01 - Tue Jul 19 13:59:22 PDT 2016
      - Analyzes a "diploid_mode 1" run of Meraculous and outputs the results

"""

#import things for rest of program
import sys
import argparse
import os
import re
import subprocess
import shutil

#https://py-gfm.readthedocs.io/en/latest/
import markdown
from mdx_gfm import GithubFlavoredMarkdownExtension

#multiprocessing stuff
from functools import partial
from multiprocessing import cpu_count
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))

def parse_arguments():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", "--censor",
                        type=str,
                        nargs='+',
                        help="""list of strings to censor for sending docs to
                        non-collaborators""")
    parser.add_argument("-q", "--quiet",
                        action='store_true')

    options = parser.parse_args()
    return options

class meraculousRunAnalyzer:
    def __init__(self, directory, censor, quiet=False):
        """In a future edit, change this method to auto-detect whether the
        run is diploid mode 1, 2, or 3
        """
        self.home = os.path.abspath(directory.strip())
        self.homeBasename = os.path.basename(directory.strip())
        self.parent = os.path.dirname(self.home.strip())

        #setup the bit to cleanly censor the files
        self.home.strip().rsplit("/")
        censor = [x for x in self.parent.strip().rsplit("/") if x]
        self.censor = {strip: "" for strip in censor}
        self.rep = dict((re.escape(k), v) for k, v in self.censor.items())
        self.pattern = re.compile("|".join(self.rep.keys()))

        #filename prefix, sanitized with censor method
        self.outname=self.str_ripper(os.path.basename(self.home))
        self.merReportsDir = os.path.join(self.parent, "meraculous_reports")
        self.reportDir = os.path.join(self.merReportsDir, self.outname)

        self.isMer = self.is_meraculous()
        self.quiet = quiet

        self.mer_size = None
        self.diploid_mode = None
        self.genome_size = None
        self.min_depth_cutoff = None
        if self.isMer:
            self.get_Params()

    def is_meraculous(self):
        os.chdir(self.home)
        subdirs = [dirs.strip() for dirs in os.listdir(self.home) if os.path.isdir(dirs.strip())]
        print(subdirs)
        os.chdir(self.parent)
        return ("log" in subdirs) and ("meraculous_import" in subdirs)

    def call_sys(self, callString, logfile):
        print("\ncalling command:", file=logfile)
        print(self.str_ripper(callString), file=logfile)
        p = subprocess.run(callString, shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True)
        output = self.str_ripper(str(p.stdout))
        err = self.str_ripper(str(p.stderr))
        print(output, file=logfile)
        print(err, file=logfile)
        if not self.quiet:
            print(output)
            print(err)
        return (output, err)

    def make_HTML(self, logfile):
        html_filepath = os.path.join(
            self.merReportsDir,
            "{0}_report.html".format(self.outname))
        print("making html at {}".format(html_filepath))
        callString = "python -m markdown -x markdown.extensions.fenced_code -x pymdownx.b64 {0} > {1}".format(
            logfile,
            html_filepath)
        p = subprocess.run(callString, shell=True)

    def get_Params(self):
        logFile = os.path.join(self.home, "log/meraculous.log")
        grep_list=["mer_size", "diploid_mode",
                   "genome_size", "min_depth_cutoff"]
        done = {x:False for x in grep_list}

        with open(logFile, "r") as f:
            for line in f:
                for query in grep_list:
                    if re.search(query, line):
                        value=line.strip().split()[1]
                        print("query: {}, value: {}".format(query, value))
                        setattr(self, query, value)
                        done[query] = True
                if False not in done.values():
                    break

    def str_ripper(self, text):
        """Got this code from here:
        http://stackoverflow.com/questions/6116978/python-replace-multiple-strings

        This method takes a set of strings, A, and removes all whole
        elements of set A from string B.

        Input: text string to strip based on instance attribute self.censor
        Output: a stripped (censored) text string
        """
        return self.pattern.sub(lambda m: self.rep[re.escape(m.group(0))], text)

    def generate_report(self):
        """This method does most of the work for this class. It generates an
        individual report for the instance.
        """

        print("making plot for {}".format(self.outname))
        # make the markdown file to log to while doing analyses
        logfile_filepath = os.path.join(
            self.reportDir,
            "{0}_report.md".format(self.outname))

        #check if seld.reportDir exists
        # if not, make it
        if not os.path.exists(self.reportDir):
            print("making the directory: {}".format(self.reportDir))
            os.makedirs(self.reportDir)

        logfile = open(logfile_filepath, 'w')

        # make a title in the markdown file
        print("# {0}".format(self.outname), file=logfile)
        print("## Meraculous assembly QC analysis", file=logfile)

        # 0. Run params
        print("### 0. Run params", file=logfile)
        print("", file=logfile)
        print("- mer_size: `{}`".format(self.mer_size), file=logfile)
        print("- diploid_mode: `{}`".format(self.diploid_mode), file=logfile)
        print("- genome_size: `{}`".format(self.genome_size), file=logfile)
        print("- min_depth_cutoff: `{}`".format(self.min_depth_cutoff), file=logfile)
        print("", file=logfile)

        # 1. get the mercount and kha plot
        print("### 1. `meraculous_mercount` output", file = logfile)
        print("", file=logfile)
        print("#### `mercount.png`", file = logfile)
        print("", file=logfile)
        mercountFromPath = os.path.join(self.homeBasename, "meraculous_mercount/mercount.png")
        mercountToPath = os.path.join(self.reportDir, "mercount.png")
        mercountHTML = os.path.join(
            os.path.basename(os.path.split(mercountToPath)[0]),
            os.path.basename(mercountToPath))
        shutil.copyfile(mercountFromPath, mercountToPath)
        print("![mercount]({0})".format(mercountHTML),
            file=logfile)
        print("", file=logfile)
        print("#### `kha.png`", file = logfile)
        print("", file=logfile)
        khaFromPath = os.path.join(self.homeBasename, "meraculous_mercount/kha.png")
        khaToPath = os.path.join(self.reportDir, "kha.png")
        khaHTML = os.path.join(
            os.path.basename(os.path.split(khaToPath)[0]),
            os.path.basename(khaToPath))
        shutil.copyfile(khaFromPath, khaToPath)
        print("![kha]({0})".format(khaHTML),
            file=logfile)
        print("", file=logfile)

        # 2. run fasta_stats on UUtig.fa
        # this is the result of the first contigs produced by meraculous_contigs
        print("### 2. `fasta_stats UUtigs.fa`", file=logfile)
        print("", file=logfile)
        print("""This is the output of `meraculous_contigs`. In a diploid genome,
              the total length of all contigs at this stage should be
              larger than the expected genome size because of the
              presence of haplotype variant UUtigs. (Meraculous
              manual)""", file=logfile)
        print("",file=logfile)
        print("```", file=logfile)
        callString="fasta_stats {0}".format(
            "{}/meraculous_contigs/UUtigs.fa".format(self.home))
        self.call_sys(callString, logfile)
        print("```", file=logfile)
        print("",file=logfile)

        if int(self.diploid_mode) == 1:
            #2.5. Check out the output of meraculous_bubble
            print("### 2.5. Diploid mode plot", file=logfile)
            print("", file=logfile)
            print("""- If there are two peaks, and the "half-depth" peak is much larger
                  than the one at full depth, that's a scenario that's
                  best handled by diploid_mode 2.  If it's the
                  opposite, make sure the parameter
                  'bubble_depth_threshold' is being set correctly (if
                  auto-detected, the result is written into
                  `meraculous_bubble/haplotigs.dmin.err`) - E. Goltsman.""",
                  file=logfile)
            print("""- For diploid assemblies, you should examine the file
                  haplotigs.depth.hist.png and verify that there are
                  two distinct peaks, one at roughly half depth of the
                  other. At this point you may want to check your
                  bubble_depth_threshold parameter and adjust it to a
                  value corresponding to the local minimum between the
                  two peaks (_ if you had originally set it to 0
                  Meraculous auto-detects this threshold_). - Meraculous Manual""",
                  file=logfile)
            print("", file=logfile)
            bblFromPath = os.path.join(self.homeBasename,
                                       "meraculous_bubble/haplotigs.depth.hist.png")
            bblToPath = os.path.join(self.reportDir, "haplotigs.depth.hist.png")
            bblHTML = os.path.join(
                os.path.basename(os.path.split(bblToPath)[0]),
                os.path.basename(bblToPath))
            shutil.copyfile(bblFromPath, bblToPath)
            print("![bbl]({0})".format(bblHTML),
                   file=logfile)
            print("", file=logfile)
            bbl_detect=os.path.join(self.home, "meraculous_bubble/haplotigs.dmin.err")
            if os.path.exists(bbl_detect):
                #in a later version, just look at the meraculous log file
                print("", file=logfile)
                print("#### `bubble_depth_threshold` auto-detected", file=logfile)
                print("", file=logfile)
                print("```", file=logfile)
                callString="cat {0}".format(bbl_detect)
                self.call_sys(callString, logfile)
                print("```", file=logfile)
                print("", file=logfile)
            else:
                print("", file=logfile)
                print("#### `bubble_depth_threshold` was user-input", file=logfile)
                print("", file=logfile)
            print("#### `fasta_stats meraculous_bubble/haplotigs.fa", file=logfile)
            print("", file=logfile)
            print("""- These are the assembly stats after running
                  `meraculous_bubble`.""", file=logfile)
            print("", file=logfile)
            print("```", file=logfile)
            callString="fasta_stats {0}".format(
            os.path.join(self.home, "meraculous_bubble/haplotigs.fa"))
            self.call_sys(callString, logfile)
            print("```", file=logfile)
            print("", file=logfile)

        if int(self.diploid_mode) == 1:
            print("### 2.6. Note on `meraculous_ono` with `diploid_mode=1`", file=logfile)
            print("", file=logfile)
            print("""- [ Note: If running in diploid_mode 1, the scaffolding is initially
                performed using combined linkage info from alternative
                variant diplotigs, i.e., both variants contribute read
                pairs to the same link as if they were one and the
                same contig. Then, using haplotype-specific read
                mapping info, phased variant paths are determined and
                the scaffold content is corrected in a haplotype
                consistent manner. As a result, one variant path
                (typically one with the higher overall depth) will be
                preserved in a multi-contig scaffold while the
                individual alternative variants are represented as
                unlinked, singleton scaffolds. A list of these
                singleton alternative variant scaffolds is also
                saved. For more on diploid-aware assembly, see section
                ‘Diploid assembly’. ] - _Meraculous manual_""", file=logfile)

        print("### 3. `meraculous_gap_closure`", file=logfile)
        print("", file=logfile)
        print("#### 3.1 `final.scaffolds.fa`", file=logfile)
        print("", file=logfile)
        print("```", file=logfile)
        callString="fasta_stats {0}".format(
            os.path.join(self.home, "meraculous_gap_closure/final.scaffolds.fa"))
        self.call_sys(callString, logfile)
        print("```", file=logfile)
        print("", file=logfile)

        if int(self.diploid_mode) == 1:
            print("#### 3.1. Note on haplotype final scaffolds", file=logfile)
            print("", file=logfile)
            print("""- [Note: when running in diploid_mode 1, if a gap
            represents a polymorphic region that had been actively
            removed earlier (i.e., a half-depth isotig), Meraculous
            will attempt to walk across it using reads from the more
            abundant allele]. - _Meraculous manual_""", file=logfile)
            print("""- When running in diploid_mode 1, a “single-haplotype” sequence
            file (final.scaffolds.single-haplotype.fa) where the
            alternative variant scaffolds have been removed is also
            created. - _Meraculous manual_""", file=logfile)
            print("", file=logfile)
            print("```", file=logfile)
            callString="fasta_stats {0}".format(
                os.path.join(self.home, "meraculous_gap_closure/final.scaffolds.single-haplotype.fa"))
            self.call_sys(callString, logfile)
            print("```", file=logfile)
            print("", file=logfile)

        print("### 4. `meraculous_final_results`", file=logfile)
        print("""- The final output of meraculous""", file=logfile)
        print("", file=logfile)
        print("```", file=logfile)
        callString="fasta_stats {0}".format(
            os.path.join(self.home, "meraculous_final_results/final.scaffolds.fa"))
        self.call_sys(callString, logfile)
        print("```", file=logfile)
        print("", file=logfile)

        # this block handles the printing the meraculous log to the new output
        # directory
        print("### 5. `meraculous_log`", file=logfile)
        print("""- The final output of meraculous""", file=logfile)
        print("", file=logfile)
        print("```", file=logfile)
        callString="cat {0}".format(
            os.path.join(self.home, "log/meraculous.log"))
        output, err = self.call_sys(callString, logfile)
        print("```", file=logfile)
        print("", file=logfile)
        new_mer_log_path = os.path.join(self.reportDir, "meraculous.log")
        with open (new_mer_log_path, "w") as f:
            print(output, file=f)

        logfile.close()
        self.make_HTML(logfile_filepath)

def run_merReport_dummy(instance):
    """This function is a helper method for merReport.generate_report().
    Specifically, it just calls the generate_report()
    function for each instance of class PreqcAnalysis to get around the
    limitations of the multiprocessing module."""
    instance.generate_report()


def determine_pool_size(job_vector):
    """This function determines how large of a pool to make based on the
    system resources currently available and how many jobs there are
    to complete.
    """

    available_threads = cpu_count()
    total_jobs = len(job_vector)
    threads_to_pass = total_jobs
    if total_jobs >= available_threads:
        threads_to_pass = available_threads
    if threads_to_pass > 90:
        threads_to_pass = 90
    print("There are {} threads available.\nThere are {} jobs:\nMaking a pool with {} threads".format(
        available_threads, total_jobs, threads_to_pass), file = sys.stderr)
    return threads_to_pass

def available_threads():
    """This is customized to always leave 6 threads on a server that has 96
    total threads, but this process takes such little time that it doesn't
    really matter.
    """
    threads = cpu_count()
    if threads >= 90:
        return 90
    else:
        return threads

def main():
    """1. Parse args
    2. Figure out which directories are actually meraculous run directories
    3. Make an instance for each directory and generate a report
    """
    home = os.getcwd()
    options = parse_arguments()
    print(options)
    print()

    #get a list of all directories in the cwd
    dirs = [direc for direc in os.listdir(home) if os.path.isdir(direc)]
    #instantiate all of the classes that we will be using in parallel processing
    instances = []
    for each in dirs:
        thisInstance = meraculousRunAnalyzer(each, options.censor, options.quiet)
        #only process the instances that are meraculous directories
        if thisInstance.isMer:
            instances.append(thisInstance)

    if len(instances) == 0:
        print("There are no meraculous folders in this directory. Exiting")
    elif len(instances) > 0:
        # run the program for each instance
        #process each file sequentially using max number of threads
        #determine the pool size to work with the unique sample names
        pool_size = determine_pool_size(instances)
        pool = ThreadPool(pool_size)
        results = pool.map(run_merReport_dummy, instances)
        pool.close()
        pool.join()

if __name__ == "__main__":
    sys.exit(main())

