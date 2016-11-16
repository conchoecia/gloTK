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

"""@author Darrin Schultz
This class tests the classes and methods for glotk_project.py
"""

import unittest
import argparse
import os
import shutil

#call the assembly with the shell
import subprocess

class assembly_test_case(unittest.TestCase):
    """Tests a swept assembly"""
    @classmethod
    def setUpClass(self):
        self.testRunDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "phix174Test")
        self.gloTKDir = os.path.join(self.testRunDir, "gloTKProjectTest")
        self.badConfigDir = os.path.join(self.testRunDir, "ksweep_test.config")
        self.goodConfigDir = os.path.join(self.testRunDir, "phix174_glob.config")


    def test_config_files(self):
        """This test checks for the exception that is raised if a symlink
        creation attempt encounters a symlink that already exists. This
        situation is caused by having the same files on two lines in the config
        file. This test also tests the converse, that no error is raised for a
        correct config file.
        """ 
        for config in [self.badConfigDir, self.goodConfigDir]:
            os.chdir(self.testRunDir)
            #make sure symlink for reads exists
            if not os.path.exists(self.gloTKDir):
                os.makedirs(self.gloTKDir)
            os.chdir(self.gloTKDir)
            callString = ["glotk-project",
                          "-i", config,
                          "-g", "pleurobrachia",
                          "-s", "bachei"]
            p = subprocess.run(callString, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)
            err = str(p.stderr)
            if config == self.badConfigDir:
                #Make sure the error is raised
                self.assertTrue("ValueError" in err)
            elif config == self.goodConfigDir:
                #check that positive test doesn't raise exception
                self.assertFalse("ValueError" in err)
            # #cleanup
            os.chdir(self.testRunDir)
            shutil.rmtree(self.gloTKDir)

    def test_files_exist(self):
        """This test simply checks to see if the requisite files are present
        after initializing a gloTK project directory.

        project_dir/
        |--glotk_info/
        |  |  default_config.config
        |  |  sample_metadata.config
        |  |  read_metadata.config
        |  |  assemblynumber_to_runname.txt
        |  |
        |  |--activity_log/
        |  |     as000.log
        |  |     as001.log
        |  |     etcetera.log
        |  |
        |  |--read_configs/
        |        reads0.yaml #this is the initial yaml file for imported reads
        |        reads1.yaml #reads generated from Seqprep, trimmomatic, et cetera
        |
        |--glotk_reads
           |  reads.log
           |
           |--reads0/
                 <forward_reads_symlink>.fq.gz
                 <reverse_reads_symlink>.fq.gz
        """
        os.chdir(self.testRunDir)
        #make sure symlink for reads exists
        if not os.path.exists(self.gloTKDir):
            os.makedirs(self.gloTKDir)
        os.chdir(self.gloTKDir)
        callString = ["glotk-project",
                      "-i", self.goodConfigDir,
                      "-g", "pleurobrachia",
                      "-s", "bachei"]
        p = subprocess.run(callString, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True)
        fileList=["gloTK_info/input_config.yaml",
                  "gloTK_info/project_init.config",
                  "gloTK_info/read_configs/reads0.yaml",
                  "gloTK_reads/reads0/SRR353630_2500_1.fastq.gz",
                  "gloTK_reads/reads0/SRR353630_2500_2.fastq.gz",
                  "gloTK_reads/reads0/SRR353630_2500_1_2.fastq.gz",
                  "gloTK_reads/reads0/SRR353630_2500_2_2.fastq.gz"]
        #verify that all of the above files exist in the correct structure
        for each in [os.path.join(self.gloTKDir, x) for x in fileList]:
            self.assertTrue(os.path.exists(each))

    def tearDown(self):
        #delete the test files once done
        os.chdir(self.testRunDir)
        if os.path.exists(self.gloTKDir):
            shutil.rmtree(self.gloTKDir)


if __name__ == '__main__':
    unittest.main()
