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

"""

title: MerParse.py
authr: darrin schultz
vrsin: 0.01

v0.01 - Tue Aug  9 18:28:25 PDT 2016
      - Started file

This program will do a few things:
1. It will read in a meraculous config file and parse the input.
2. It will write new config files from the input and perform a multi-assembly
   parameter sweep. Simply put, it runs many assemblies through a user-defined
   parameter range so that the user can find the optimum assembly.
"""
from .libseq import LibSeq
from collections import UserDict
from numbers import Number
from time import strftime as tfmt
from string import ascii_letters, digits
import os

class MerParse:
    """This class:

    1. reads in an input file and extrapolates all of the possible
       parameters for Meraculous
    2. Choose a parameter to sweep on and calculates the sweep intervals
    3. Output all of the config files to <$PWD>/mer_configs
    4. Return a list of filepaths for the main method to use for assembly and
        report generation.

    Useage example:
    myMerParser = MerParse(<params>)
    myPaths = myMerParser.sweeper_output()

    Now the myPaths object contains a dict of run names and absolute paths for
        the config files.
    """
    def __init__(self, inputFile, sweep, sStart, sStop, sInterval, lnProcs,
                 asPrefix = "as", asSI = 0, genus = None, species = None):
        # --------------- Config File Parameters -------------------------------
        self.inputFile = inputFile
        self.params = {"lib_seq": [],
                       "genome_size": -0.1,
                       "mer_size": -1,
                       "min_depth_cutoff": 0,
                       "num_prefix_blocks": -1,
                       "diploid_mode": 0,
                       "use_cluster": 0,
                       "no_read_validation": 0,
                       "fallback_on_est_insert_size": 0,
                       "gap_close_aggressive": 0,
                       "gap_close_rpt_depth_ratio": 2.0,
                       "local_num_procs": 1,
                       "local_max_retries": 0}
        self.diploid_mode = {"bubble_depth_threshold": 0,
                             "strict_haplotypes": 1}
        # makes two lists that contain which strings refer to things that should
        #  be converted to ints vs floats
        self.to_int = [x for x in self.params if isinstance(self.params[x], int)] + [
            y for y in self.diploid_mode if isinstance(self.diploid_mode[y], int)]
        self.to_float = [x for x in self.params if isinstance(self.params[x], float)] + [
            y for y in self.diploid_mode if isinstance(self.diploid_mode[y], float)]

        # for debugging purposes, makes a list of values that were assigned a
        #  negative value upon initialization, showing that they were not
        #  specified in the config file and need to be changed. Raise an error
        #  in the read_config file.
        # self.config_specified are things that must be specified in the config
        #  file, so if they are negative we know that the user made a mistake.
        # self.config_unspecified are things that have a default value already,
        #  so if they are negative then the user mis-entered something in the
        # config file.
        self.config_specified = [x for x in [u for u in self.params if isinstance(self.params[u], Number)] if self.params[x] < 0]
        self.config_unspecified = [x for x in [u for u in self.params if isinstance(self.params[u], Number)] if self.params[x] >= 0] + [
            y for y in [v for v in self.diploid_mode if isinstance(self.diploid_mode[v], Number)] if self.diploid_mode[y] >= 0]

        # -------------------------Sweep Parameters-----------------------------
        self.sweep = sweep
        self.sweep_support_int = ["mer_size", "bubble_depth_threshold"]
        self.sStart = sStart
        self.sStop = sStop
        self.sInterval = sInterval

            #make sure that they're supported sweep types
        if self.sweep not in self.sweep_support_int:
            raise AttributeError("""{0} is not a currently supported sweep-able
            parameter. Please choose one of the following parameters or update
            the software at https://github.com/cypridina/gloTK. {1}""".format(
                self.sweep, self.sweep_support))

            #if they're supposed to be ints, convert them to ints first
        if self.sweep in self.sweep_support_int:
            for each in ["sStart", "sStop", "sInterval"]:
                setattr(self, each, int(getattr(self, each)))
                attr=getattr(self, each)
                if attr % 1 != 0:
                    raise AttributeError("""One or more of your sweep parameters
                    contains float numbers {0}, but you chose to sweep on a parameter
                    that requires ints {1}. Please try again.""".format(each, self.sweep))

        if self.sStart >= self.sStop:
            raise AttributeError("""The input sweep start value is greater than
            or equal to the sweep stop value. {0}>={1}. Please make sweep start
            less than sweep stop.""".format(self.sStart, self.sStop))
        if ((self.sStop - self.sStart) % self.sInterval) != 0:
            raise AttributeError("""Your sweep interval does not end on your
            sweep stop, and is off by {0}. Try increasing or decreasing your
            sweep stop by that much.""".format(self.sInterval - ((self.sStop - self.sStart) % self.sInterval)))
                #------------------------------Sweep List-------------------------------
        # Generate list for the sweep parameters. Current implementation only
        #   sweeps through ints, so no special processing is needed.
        # If sweep is "mer_size", then make sure that none of the values are
        #   even.
        if self.sweep in self.sweep_support_int:
            sweep_list = [x for x in range(self.sStart, self.sStop + self.sInterval, self.sInterval)]
            if self.sweep == "mer_size":
                evens = [x for x in sweep_list if x % 2 == 0]
                if evens:
                    raise AttributeError("""Your sweep interval for mer_size
                    resulted in the following even numbers: {0}. Please adjust
                    your sweep start, stop, and interval to only contain
                    positive odd integers""".format(evens))
            self.sList = sweep_list

        # -------------------------Name Parameters------------------------------
        self.as_a = asPrefix
        self.as_i = int(asSI)
        #http://stackoverflow.com/questions/311627/
        self.as_d = tfmt("%Y%m%d")
        self.as_z = "ME"
        self.find_illegal_characters(genus, species)
        self.as_g = genus if genus else None
        self.as_s = species if species else None
        # local number of processors (number of processors for one run)
        self.lnProcs = lnProcs

        #-----------------------Directory Parameters----------------------------
        self.cwd = os.getcwd()
        #-----------------------Read the config file----------------------------
        self.read_config()


    def find_illegal_characters(self, genus, species):
        #only allow the user to input ASCII letters and digits (no punctuation)
        # in the genus and species names to avoid filename errors
        illegal = [x for x in
                   str(genus).replace("None","") + str(species).replace("None","")
                   if x not in ascii_letters + digits]
        if illegal:
            raise IOError("""<genus> or <species/sample> contains the
                illegal characters: {}""".format(illegal))

    def name_gen(self, as_number, k):
        """See the sweeper_output() method doctstring for format details."""
        #how to add prefixes to numbers
        #http://stackoverflow.com/a/135157
        #         ai_d_z_g_s_k

        parts =  {"1ai": "{0}{1:03d}_".format(self.as_a, as_number),
                  "2d" : "{0}_".format(self.as_d),
                  "3z" : "ME_",
                  "4g" : "{}_".format(self.as_g) if self.as_g else "",
                  "5s" : "{}_".format(self.as_s) if self.as_s else "",
                  "6k" : "k{}".format(k)}
        return "".join([parts.get(x) for x in sorted(parts)])

    def sweeper_output(self):
        """This method outputs all of the config files in the configs directory
        within the current working directory. It names the config files based on
        the following format:

        a = assembly prefix <default = as>
        i = assembly start index <default = 000>
        d = date in format YYYYMMDD <default=today>
        z = assembler type - ME for Meraculous in this case
        g = first four letters of genus. Truncates if species name is too long.
        s = first four letters of species. Truncates if genus name is too long.
        k = kmer size

        ai_d_z_g_s_k

        For example, if I enter the following parameters for a k=21 run on
        20160811:
          - g = "Malacosteus"
          - s = "niger"

        The assembly directory will be saved in the current working
        directory as:
          - configs/as000_20160811_ME_mala_nige_k21

        Steps for this method:
        1. Check if configs is a directory that exists, if not, make it.
        2. For loop
          - Assign assembly numbers to sweep param
          - Output config file for each sweep param
        3. return dict of {"<run string>": "<config abs path>"}
        """
        # 1. check if configs directory exists, if not, make it
        config_dir = os.path.join(os.getcwd(), "configs")
        if not os.path.exists(config_dir):
            os.makedirs(config_dir)

        # 2. Assign assembly numbers to sweep param
        asNum_sweep_dict = {asNum:self.sList[asNum - self.as_i]
                            for asNum in range(self.as_i, self.as_i + len(self.sList))}
        if self.sweep == "mer_size":
            asName_asSweep_dict = {
                self.name_gen(asNum, asNum_sweep_dict[asNum]): asNum_sweep_dict[asNum]
                                      for asNum in asNum_sweep_dict}
        else:
            asName_asSweep_dict = {
                self.name_gen(asNum, self.params.get("mer_size")): asNum_sweep_dict[asNum]
                                  for asNum in asNum_sweep_dict}

        for name in asName_asSweep_dict:
            sweep = asName_asSweep_dict[name]
            with open(os.path.join(config_dir, name), "w") as f:
                self.params[self.sweep] = sweep
                for key in self.params:
                    if key == "lib_seq":
                        for each in self.params["lib_seq"]:
                            print("lib_seq {0}".format(each), file = f)
                    else:
                        print("{0} {1}".format(key,self.params[key]), file=f)

        # 3 return dict of {"<run string>": "<config abs path>"}
        return {name:os.path.join(config_dir, name)
                for name in asName_asSweep_dict}

    def assign(self, dict_name, key, value):
        """This assigns an input string to either a float or int and saves it in
        the config params as such"""
        if key in self.to_int:
            getattr(self, dict_name)[key] = int(value)
        elif key in self.to_float:
            getattr(self, dict_name)[key] = float(value)
        else:
            raise ValueError("""The param:value pair, {0}:{1}, has neither been
            identified as an int nor a float. Therefore, we are not able to save
            its value to any dictionary in the merParse class""".format(
                key, value))

    def libseq_parse(self, line):
        """This method parses a line containing a lib_seq string and adds it to
        the lib_seq key in the self.params dictionary"""

        LS = LibSeq()
        line = line.split()
        if (len(line) - 1) != 12:
            raise ValueError("""ERROR: There are {0} values in the lib_seq line when
                             there should only be 10. Check your input meraculous
                             config file for errors and make sure that it matches
                             the specification for the Meraculous manual.""".format(
                                 len(line)))
        for index in LS.values:
            LS[LS.values.get(index)] = line[index]
        return LS

    def space_error(self, line):
        raise ValueError(
            """ERROR: There are too many spaces in this line in your config
            file: {0}""".format(line))

    def read_config(self):
        with open(self.inputFile, "r") as f:
            for line in f:
                line=line.strip()
                # if the line is empty, this will prevent observing those lines
                if line:
                    if line[0] is not "#":
                        split = line.split()
                        param = split[0]
                        val = split[1]
                        if param in self.params:
                            if param == "lib_seq":
                                #add a libseq object to the libseq key in self.params
                                self.params["lib_seq"].append(self.libseq_parse(line))
                            else:
                                if len(split) == 2:
                                    self.assign("params", param, val)
                                else:
                                    #make sure there is only a param and value here,
                                    # otherwise there must be an error in the config
                                    # file on this line
                                    self.space_error(line)

                        elif param in self.diploid_mode:
                            if len(split) == 2:
                                self.assign("diploid_mode", param, val)
                            else:
                                # should only be 2 elements long
                               self.space_error(line)
                        # If the param isn't in either of those two dictionaries, it
                        # hasn't been implemented yet
                        else:
                            print("""This line is not legal input for this version of
                            glotk-mer-ksweep.py. Please remove it or consult
                            the software author via
                            https://github.com/cypridina/gloTK

                            {0}""".format(line))
        # check to make sure that the following things are not negative values,
        # signalling that they were never appropriately changed
        for tempDict in [self.params, self.diploid_mode]:
            for key in tempDict:
                val = tempDict.get(key)
                # This avoids comparasion problems with the LibSeq objects
                if key != "lib_seq":
                    if (val < 0) and (key in self.config_specified):
                        #see the __init__() method at the self.config_specified
                        # definition for a better description
                        raise ValueError(
                            """The param {0} was not specified in the config file.
                            Please specify it in the config file and re-run this
                            script""".format(key))
                    elif (val < 0) and (key in self.config_unspecified):
                        raise ValueError(
                            """The param {0} was incorrectly specified as a negative
                            number in the config file. Please use only positive
                            integers or positive floats in the config file and try
                            again.""".format(key))
        # override the number of procs based on parallelism controlled in
        #  gloTK/scripts/glotk_sweep.py
        self.params["local_num_procs"] = self.lnProcs
