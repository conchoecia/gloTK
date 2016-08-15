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
from glotk_mer_ksweep import libSeq, merParse
from time import strftime as tfmt
import re
import ast


class libSeq_test_case(unittest.TestCase):
    """Tests that libSeq works correctly, like a dictionary
    """
    def test_libSeq_construction(self):
        myLib = libSeq()
        values = ["wildcard", "name", "insertAvg", "insertSdev", "avgReadLn",
                  "hasInnieArtifact", "isRevComped", "useForContiging",
                  "scaffRound", "useForGapClosing",
                  "5p_wiggleRoom", "3p_wiggleRoom"]
        data = {x: "" for x in values}
        for key in myLib.keys():
            data.pop(key)
        self.assertEqual(data, {})

class merParse_test_case(unittest.TestCase):
    """Tests that the merParse class works correctly"""
    def setUp(self):
        self.configPath = "ksweep_test.config"
        self.sweep = "mer_size"
        self.sStart = 21
        self.sStop = 75
        self.sInterval = 2
        self.asPrefix = "ab"
        self.asSI = 5
        self.genus = "Malacosteus"
        self.species = "niger"
        self.myParse = merParse(self.configPath, self.sweep,
                                self.sStart, self.sStop,
                                self.sInterval)

        self.correctly_parsed_ls1 = {"wildcard": "/mydir/fwd_data*.1.fastq.gz,/mydir/rev_data*.2.fastq.gz",
                            "name": "SP2013",
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

        self.correctly_parsed_ls2 = {"wildcard": "/mydir1*.fastq.gz,/mydir2*.fastq.gz",
                            "name": "SPAGET",
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

    def test_libseq_parse(self):
        """This tests if the libseq_parse() method correctly parses "lib_seq"
        lines in the config files"""

        line = "lib_seq /mydir/fwd_data*.1.fastq.gz,/mydir/rev_data*.2.fastq.gz SP2013 720 100 100 0 0 1 1 1 0 0"
        parsed = self.myParse.libseq_parse(line)
        self.assertEqual(parsed, self.correctly_parsed_ls1)

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
        self.assertEqual(self.myParse.params, parsed_params)
        self.assertEqual(self.myParse.diploid_mode, parsed_diploid_mode)

    def test_sweepStart_less(self):
        """For this program, the sweep start must be less than sweep stop,
        so it throws an exception when that rule is violated

        Tests this line of code in merParse class:
          - if sStart >= sStop:
        """
        with self.assertRaises(AttributeError):
            self.myParse = merParse(self.configPath, self.sweep,
                                self.sStop, self.sStart, self.sInterval)

    def test_sweep_supported(self):
        """Checks to make sure that an error is raised if something that isn't
        one of the supported sweep types is input as self.sweep in merParse
        class.

        Testing block that starts with:
          - if self.sweep not in self.sweep_support_int:
        """
        #shouldn't work with an unknown string
        with self.assertRaises(AttributeError):
            self.myParse = merParse(self.configPath, "noodles",
                                self.sStop, self.sStart, self.sInterval)
        #and it shouldn't work with a number
        with self.assertRaises(AttributeError):
            self.myParse = merParse(self.configPath, 5,
                                self.sStop, self.sStart, self.sInterval)

    def test_sweep_mer_size_odd(self):
        """tests to make sure that all off the mer_sizes for sweep are odd
        positive integers. This also handles the case that the start and stop
        are entered as even numbers since the error prints out which numbers are
        even.

        Tests this block of code:
          - if evens:
        """
        with self.assertRaises(AttributeError):
            self.myParse = merParse(self.configPath, 5, 21, 27, 3)

    ############################################################################
    #                          myParse.name_gen() Tests
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
        myParseShouldEqual = "as000_{}_ME_k21".format(as_d)
        equals = self.myParse.name_gen(0, 21)
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
        myParseShouldEqual = "ab015_{}_ME_Malacosteus_k{}".format(as_d, k)
        self.myParse2 = merParse(self.configPath, self.sweep,
                                self.sStart, self.sStop,
                                self.sInterval, asPrefix = "ab",
                                genus = self.genus, asSI=5)
        self.assertEqual(myParseShouldEqual, self.myParse2.name_gen(15, k))

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
        myParseShouldEqual = "xx856_{}_ME_niger_k{}".format(as_d, k)
        self.myParse2 = merParse(self.configPath, self.sweep,
                                self.sStart, self.sStop,
                                self.sInterval, asPrefix = "xx",
                                species = self.species, asSI=5)
        self.assertEqual(myParseShouldEqual, self.myParse2.name_gen(856, k))

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
            self.myParse3 = merParse(self.configPath, self.sweep,
                                self.sStart, self.sStop,
                                self.sInterval, asPrefix = "xx",
                                species = species,
                                genus = genus, asSI=5)


        exception = raises_cm.exception
        string_list = re.findall( "\[.+", exception.args[0])[0]
        list_list = ast.literal_eval(string_list)
        self.assertEqual(list_list, illegal_chars)
    ############################################################################
    #                      DONE myParse.name_gen() Tests DONE
    ############################################################################

    ############################################################################
    #                       myParse.sweeper_output() Tests
    ############################################################################
    @unittest.skip("only run this if you want to see output files in your terminal")
    def test_sweeper_output(self):
        myParseSweeper = merParse(self.configPath, self.sweep,
                                self.sStart, self.sStop,
                                self.sInterval, asPrefix = "xx",
                                species = self.species,
                                genus = self.genus, asSI=5)
        myParseSweeper.sweeper_output()


if __name__ == '__main__':
    unittest.main()
