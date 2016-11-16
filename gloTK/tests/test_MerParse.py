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
This class tests the classes and methods for glotk-mer-ksweep.py
"""

import unittest
from gloTK import MerParse, LibSeq, ConfigParse
from time import strftime as tfmt
import os
import re
import ast

class merParse_test_case(unittest.TestCase):
    """Tests that the merParse class works correctly"""
    def setUp(self):
        self.configPath = os.path.join(os.path.abspath(os.path.dirname(__file__)),"phix174Test/ksweep_test.config")
        self.testParser = ConfigParse(self.configPath)
        self.pwd = os.path.dirname(os.path.abspath(__file__))

        self.sweep = "mer_size"
        self.sList = ['21', '23', '25', '27', '29',
                      '31', '33', '35', '37', '39',
                      '41', '43', '45', '47', '49',
                      '51', '53', '55', '57', '59',
                      '61', '63', '65', '67', '69',
                      '71', '73', '75']
        self.asPrefix = "ab"
        self.asSI = 5
        self.genus = "Malacosteus"
        self.species = "niger"
        self.lnProcs = 27
        self.myParse = MerParse(self.configPath, self.sweep,
                                self.sList,
                                self.lnProcs)

        self.correctly_parsed_ls1 = {"wildcard": "reads/SRR353630_2500_1*.fastq.gz,reads/SRR353630_2500_2*.fastq.gz",
                            "name": "SP2013",
                            "globs":[os.path.join(self.pwd, "phix174Test/reads/SRR353630_2500_1*.fastq.gz"),
                                     os.path.join(self.pwd, "phix174Test/reads/SRR353630_2500_2*.fastq.gz")],
                            "pairs":[(os.path.join(self.pwd, 'phix174Test/reads/SRR353630_2500_1.fastq.gz'),
                                     os.path.join(self.pwd, 'phix174Test/reads/SRR353630_2500_2.fastq.gz')),

                                     (os.path.join(self.pwd, 'phix174Test/reads/SRR353630_2500_1_2.fastq.gz'),
                                     os.path.join(self.pwd, 'phix174Test/reads/SRR353630_2500_2_2.fastq.gz'))], 
                            "insertAvg": "720",
                            "insertSdev": "100" ,
                            "avgReadLn": "100",
                            "hasInnieArtifact": "0",
                            "isRevComped": "0",
                            "useForContiging": "1",
                            "scaffRound": "1",
                            "useForGapClosing": "1",
                            "5p_wiggleRoom": "0",
                            "3p_wiggleRoom": "0"}


        self.correctly_parsed_ls2 = {"wildcard": "reads/SRR353630_2500_1*.fastq.gz,reads/SRR353630_2500_2*.fastq.gz",
                            "name": "SPAGET",
                            "globs":[os.path.join(self.pwd, "phix174Test/reads/SRR353630_2500_1*.fastq.gz"),
                                     os.path.join(self.pwd, "phix174Test/reads/SRR353630_2500_2*.fastq.gz")],
                            "pairs":[(os.path.join(self.pwd, 'phix174Test/reads/SRR353630_2500_1.fastq.gz'),
                                     os.path.join(self.pwd, 'phix174Test/reads/SRR353630_2500_2.fastq.gz')),

                                     (os.path.join(self.pwd, 'phix174Test/reads/SRR353630_2500_1_2.fastq.gz'),
                                     os.path.join(self.pwd, 'phix174Test/reads/SRR353630_2500_2_2.fastq.gz'))], 

                            "insertAvg": "555",
                            "insertSdev": "400" ,
                            "avgReadLn": "150",
                            "hasInnieArtifact": "1",
                            "isRevComped": "1",
                            "useForContiging": "0",
                            "scaffRound": "0",
                            "useForGapClosing": "0",
                            "5p_wiggleRoom": "1",
                            "3p_wiggleRoom": "1"}

    def test_read_config(self):
        parsed_params = {"lib_seq": [self.correctly_parsed_ls1, self.correctly_parsed_ls2],
                       "genome_size": 0.328,
                       "mer_size": 71,
                       "min_depth_cutoff": 0,
                       "num_prefix_blocks": 1,
                       "diploid_mode": 1,
                       "use_cluster": 0,
                       "no_read_validation": 0,
                       "fallback_on_est_insert_size": 0,
                       "gap_close_aggressive": 0,
                       "gap_close_rpt_depth_ratio": 2.0,
                       "local_num_procs": 27,
                       "local_max_retries": 0}
        parsed_diploid_mode = {"bubble_depth_threshold": 0,
                             "strict_haplotypes": 1}
        self.maxDiff = None
        self.assertEqual(self.testParser.params, parsed_params)
        self.assertEqual(self.testParser.diploid_mode, parsed_diploid_mode)

    def test_sweep_supported(self):
        """Checks to make sure that an error is raised if something that isn't
        one of the supported sweep types is input as self.sweep in merParse
        class.

        Testing block that starts with:
          - if self.sweep not in self.sweep_support_int:
        """
        #shouldn't work with an unknown string
        with self.assertRaises(AttributeError):
            self.myParse = MerParse(self.configPath, "noodles",
                                self.sList,
                                self.lnProcs)
        #and it shouldn't work with a number
        with self.assertRaises(AttributeError):
            self.myParse = MerParse(self.configPath, 5,
                                self.sList,
                                self.lnProcs)

    def test_sweep_mer_size_odd(self):
        """tests to make sure that all off the mer_sizes for sweep are odd
        positive integers. This also handles the case that the start and stop
        are entered as even numbers since the error prints out which numbers are
        even.

        Tests this block of code:
          - if evens:
        """
        with self.assertRaises(AttributeError):
            self.myParse = MerParse(self.configPath, "mer_parse", ['21', '22'],
                                    self.lnProcs)

    ############################################################################
    #                          MerParse.name_gen() Tests
    ############################################################################
    def test_name_gen_basic(self):
        """This tests that the myParse.name_gen() method works in instances
        where the user neither defines the input genus nor species. The output
        should just have an omitted genus and species with no extra underscores
        """
        # myParse example
        #   - init params were:
        #     - asPrefix = "ab"
        #     - genus = None
        #     - species = None
        #     - assembly_prefix = as (default)
        #     - assembly_startI = 0 (default)
        #     - k = 21
        # should be "as000_YYYYMMDD_ME_k21
        as_d = tfmt("%Y%m%d")
        myParseShouldEqual = "as000_{}_ME_k21_d2".format(as_d)
        equals = self.myParse.name_gen(0, 21, 2)
        self.assertEqual(myParseShouldEqual, equals)

    def test_name_gen_noSpecies(self):
        """This tests that the myParse.name_gen() method works when the
        species is not defined. The output should just have an omitted
        genus with no extra underscores.
        """
        # myParse2 example
        #   - init params were:
        #     - asPrefix = "ab"
        #     - genus = Malacosteus
        #     - species = None
        #     - assembly selected i = 15
        #     - k = 33
        # should be "as015_YYYYMMDD_ME_k33
        as_d = tfmt("%Y%m%d")
        k = 33
        myParseShouldEqual = "ab015_{}_ME_Malacosteus_k{}_d0".format(as_d, k)
        self.myParse2 = MerParse(self.configPath, self.sweep,
                                self.sList, self.lnProcs,
                                asPrefix = "ab",
                                genus = self.genus, asSI=5)
        self.assertEqual(myParseShouldEqual, self.myParse2.name_gen(15, k, 0))

    def test_name_gen_noGenus(self):
        """This tests that the myParse.name_gen() method works when the
        genus is not defined. The output should just have an omitted
        specieswith no extra underscores.
        """
        # myParse3 example
        #   - init params were:
        #     - asPrefix = "xx"
        #     - genus = None
        #     - species = niger
        #     - assembly selected i = 856
        #     - k = 57
        # should be "as856_YYYYMMDD_ME_niger_k57
        as_d = tfmt("%Y%m%d")
        k = 57
        myParseShouldEqual = "xx856_{}_ME_niger_k{}_d0".format(as_d, k)
        self.myParse2 = MerParse(self.configPath, self.sweep,
                                self.sList, self.lnProcs,
                                asPrefix = "xx",
                                species = self.species, asSI=5)
        self.assertEqual(myParseShouldEqual, self.myParse2.name_gen(856, k, 0))

    def test_genus_species_illegal(self):
        """This tests that the myParse.__init__() method correctly identifies
        the illegal characters (non-ascii and non-digits) in the genus
        and species/sample parameters.

        """
        # myParse3 example
        #   - init params were:
        #     - asPrefix = "xx"
        #     - genus = None
        #     - species = niger
        #     - assembly selected i = 856
        #     - k = 57
        # should be "as856_YYYYMMDD_ME_niger_k57
        as_d = tfmt("%Y%m%d")
        k = 57
        myParseShouldEqual = "xx856_{}_ME_niger_k{}".format(as_d, k)
        species = "Yello`~!@#$%^&*()_-+={[}]|\\\"':;?/>.<,"
        genus = "CAT"
        illegal_chars = ["`", "~", "!", "@", "#", "$", "%", "^", "&", "*", "(",
                         ")", "_", "-", "+", "=", "{", "[", "}", "]", "|", "\\",
                         "\"","'", ":", ";", "?", "/", ">", ".", "<", ","]
        #http://stackoverflow.com/questions/12516881
        with self.assertRaises(OSError) as raises_cm:
            self.myParse3 = MerParse(self.configPath, self.sweep,
                                self.sList, self.lnProcs,
                                asPrefix = "xx",
                                species = species,
                                genus = genus, asSI=5)


        exception = raises_cm.exception
        string_list = re.findall( "\[.+", exception.args[0])[0]
        list_list = ast.literal_eval(string_list)
        self.assertEqual(list_list, illegal_chars)
    ############################################################################
    #                      DONE MerParse.name_gen() Tests DONE
    ############################################################################

    ############################################################################
    #                       MerParse.sweeper_output() Tests
    ############################################################################
    @unittest.skip("only run this if you want to see output files in your terminal")
    def test_sweeper_output(self):
        myParseSweeper = MerParse(self.configPath, self.sweep,
                                self.sList, asPrefix = "xx",
                                species = self.species,
                                genus = self.genus, asSI=5)
        myParseSweeper.sweeper_output()


if __name__ == '__main__':
    unittest.main()
