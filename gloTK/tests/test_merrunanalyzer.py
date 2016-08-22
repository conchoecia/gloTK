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
This class tests the methods for the MerRunAnalyzer class. Important things to
test include if the class correctly moves files and generates a readable output.
"""
from gloTK import MerRunAnalyzer

import unittest
import argparse
import os
import shutil

class merParse_test_case(unittest.TestCase):
    """Tests that the merParse class works correctly"""

    def setUp(self):
        self.testRunDir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "meraculousTestRun")
        self.outputParentDir = os.path.join(os.path.abspath(os.path.dirname(__file__)))
        self.censor = []

    def tearDown(self):
        """This directory is called after every method is run even if the method
        produces an error.

        Here, it is used to remove any report files that are generated while
        running these tests.
        """
        self.outputParentDir = os.path.join(os.path.abspath(os.path.dirname(__file__)) ,
                                            "reports")
        if os.path.exists(self.outputParentDir):
            shutil.rmtree(self.outputParentDir)


    def test_is_meraculous(self):
        """Tests that the software correctly detects if this is a meraculous
        directory. This doesn't require that anything further than the log file
        was created and meraculous_import started."""
        reporter = MerRunAnalyzer(self.testRunDir, self.outputParentDir,
                                  self.censor)
        self.assertTrue(reporter.isMer)

    def test_is_not_meraculous(self):
        """Tests that the software correctly detects if something isn't a
        meraculous directory."""
        self.notMerDir = os.path.join(os.path.abspath(os.path.dirname(__file__)))
        reporter = MerRunAnalyzer(self.notMerDir, self.outputParentDir,
                                  self.censor)
        self.assertFalse(reporter.isMer)

    def test_html_built(self):
        """Checks to verify if the HTML file and the directory containing the
        plots, etc were created.

        The directory structure for the reports is:

        reports/
        |  meraculousTestRun_report.html
        |  meraculousTestRun/
        |  |  haplotigs.depth.hist.png
        |  |  kha.png
        |  |  meraculous.log
        |  |  meraculousTestRun_report.md
        |  |  mercount.png
        """
        reporter = MerRunAnalyzer(self.testRunDir, self.outputParentDir,
                                  self.censor)
        reporter.generate_report()
        #The parent report directory
        self.outputParentDir = os.path.join(os.path.abspath(os.path.dirname(__file__)) ,
                                            "reports")
        self.assertTrue(os.path.isdir(self.outputParentDir))

        #The enclosed HTML report file
        self.htmlReport = os.path.join(self.outputParentDir, "meraculousTestRun_report.html") 
        self.assertTrue(os.path.exists(self.htmlReport))

        #The enclosed report documents directory
        self.reportDocs = os.path.join(self.outputParentDir, "meraculousTestRun")
        self.assertTrue(os.path.isdir(self.reportDocs))

        #The enclosed diploid mode PNG
        self.bubblePNG = os.path.join(self.reportDocs, "haplotigs.depth.hist.png")
        self.assertTrue(os.path.exists(self.bubblePNG))

        #the kha plot
        self.khaPNG = os.path.join(self.reportDocs, "kha.png")
        self.assertTrue(os.path.exists(self.khaPNG))

        #the meraculous log
        self.merLog = os.path.join(self.reportDocs, "meraculous.log")
        self.assertTrue(os.path.exists(self.merLog))

        #the meraculous raw .md file
        self.md = os.path.join(self.reportDocs, "meraculousTestRun_report.md")
        self.assertTrue(os.path.exists(self.md))

        #the mercount png
        self.mercountPNG = os.path.join(self.reportDocs, "mercount.png")
        self.assertTrue(os.path.exists(self.mercountPNG))

        #remove everything after d
if __name__ == '__main__':
    unittest.main()
