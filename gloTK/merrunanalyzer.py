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


class MerRunAnalyzer:
    """This class generates a report for one meraculous run. It requires a run
    directory to look for the report, as well as a directory to deposit the
    report. The generate_report() method will output an HTML file and folder
    in the report directory. It will use the meraculous run name as the name for
    the report HTML file and report directory.
    """
    def __init__(self, run_directory, output_parent_dir,
                 censor, quiet=False):
        """In a future edit, change this method to auto-detect whether the
        run is diploid mode 1, 2, or 3
        """

        self.home = os.path.abspath(run_directory.strip())
        self.reportName = os.path.basename(run_directory.strip())
        self.parent = os.path.dirname(self.home.strip())

        #setup the bit to cleanly censor the files
        self.home.strip().rsplit("/")
        censor = [x for x in self.parent.strip().rsplit("/") if x]
        self.censor = {strip: "" for strip in censor}
        self.rep = dict((re.escape(k), v) for k, v in self.censor.items())
        self.pattern = re.compile("|".join(self.rep.keys()))

        #filename prefix, sanitized with censor method
        self.outname=self._str_ripper(os.path.basename(self.home))
        self.merReportsDir = os.path.join(self.parent, "reports")
        self.reportDir = os.path.join(self.merReportsDir, self.outname)

        #make sure that the directory is a meraculous directory
        self.isMer = self._is_meraculous()
        self.quiet = quiet

        self.mer_size = None
        self.diploid_mode = None
        self.genome_size = None
        self.min_depth_cutoff = None
        if self.isMer:
            self._get_Params()

    def _is_meraculous(self):
        subdirs = [dirs.strip() for dirs in os.listdir(self.home) if os.path.isdir(os.path.join(self.home,dirs.strip()))]
        return ("log" in subdirs) and ("meraculous_import" in subdirs)

    def _call_sys(self, callString, logfile):
        print("\ncalling command:", file=logfile)
        print(self._str_ripper(callString), file=logfile)
        p = subprocess.run(callString, shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True)
        output = self._str_ripper(str(p.stdout))
        err = self._str_ripper(str(p.stderr))
        print(output, file=logfile)
        print(err, file=logfile)
        # if not self.quiet:
        #     print(output)
        #     print(err)
        return (output, err)

    def _make_HTML(self, logfile):
        html_filepath = os.path.join(
            self.merReportsDir,
            "{0}_report.html".format(self.outname))
        print("making html at {}".format(html_filepath))
        callString = "python -m markdown -x markdown.extensions.fenced_code -x pymdownx.b64 {0} > {1}".format(
            logfile,
            html_filepath)
        p = subprocess.run(callString, shell=True)

    def _get_Params(self):
        logFile = os.path.join(self.home, "log/meraculous.log")
        grep_list=["mer_size", "diploid_mode",
                   "genome_size", "min_depth_cutoff"]
        done = {x:False for x in grep_list}

        with open(logFile, "r") as f:
            for line in f:
                for query in grep_list:
                    if re.search(query, line):
                        value=line.strip().split()[1]
                        #print("query: {}, value: {}".format(query, value))
                        setattr(self, query, value)
                        done[query] = True
                if False not in done.values():
                    break

    def _str_ripper(self, text):
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

        print("making HTML report for {}".format(self.outname))
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
        mercountFromPath = os.path.join(self.home, "meraculous_mercount/mercount.png")
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
        khaFromPath = os.path.join(self.home, "meraculous_mercount/kha.png")
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
        self._call_sys(callString, logfile)
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
            bblFromPath = os.path.join(self.home,
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
                self._call_sys(callString, logfile)
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
            self._call_sys(callString, logfile)
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
        self._call_sys(callString, logfile)
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
            self._call_sys(callString, logfile)
            print("```", file=logfile)
            print("", file=logfile)

        print("### 4. `meraculous_final_results`", file=logfile)
        print("""- The final output of meraculous""", file=logfile)
        print("", file=logfile)
        print("```", file=logfile)
        callString="fasta_stats {0}".format(
            os.path.join(self.home, "meraculous_final_results/final.scaffolds.fa"))
        self._call_sys(callString, logfile)
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
        output, err = self._call_sys(callString, logfile)
        print("```", file=logfile)
        print("", file=logfile)
        new_mer_log_path = os.path.join(self.reportDir, "meraculous.log")
        with open (new_mer_log_path, "w") as f:
            print(output, file=f)

        logfile.close()
        self._make_HTML(logfile_filepath)

