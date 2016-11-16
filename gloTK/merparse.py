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
import copy
import os
import yaml

class MerParse:
    """This class:

    1. reads in an input file into ConfigParse and extrapolates all of the possible
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
    def __init__(self, inputFile, sweep=None, sList=None, lnProcs=0,
                 asPrefix = "as", asSI = 0, genus = None, species = None,
                 triplet=False):
        # ----------------- Run triplet assembly -------------------------------
        self.triplet = triplet
        # --------------- Config File Parameters -------------------------------
        self.inputFile = inputFile
        configFile = ConfigParse(inputFile)
        self.params = configFile.params
        self.diploid_mode = configFile.diploid_mode
        # override the number of procs based on parallelism controlled in
        #  gloTK/scripts/glotk_sweep.py
        self.params["local_num_procs"] = lnProcs


       # -------------------------Sweep Parameters-----------------------------
        self.sweep = sweep
        self.sweep_support_int = ["mer_size", "bubble_depth_threshold"]

            #make sure that they're supported sweep types
        if self.sweep not in self.sweep_support_int:
            raise AttributeError("""{0} is not a currently supported sweep-able
            parameter. Please choose one of the following parameters or update
            the software at https://github.com/cypridina/gloTK. {1}""".format(
                self.sweep, self.sweep_support))


        #------------------------------Sweep List-------------------------------
        # Generate list for the sweep parameters. Current implementation only
        #   sweeps through ints, so no special processing is needed.
        # If sweep is "mer_size", then make sure that none of the values are
        #   even.
        if self.sweep in self.sweep_support_int:
            if self.sweep == "mer_size":
                evens = [ x for x in sList if int(x) % 2 == 0]
                floats = [ x for x in sList if float(x) % int(x) != 0.0 ]
                errors = evens + floats
                if errors:
                    raise AttributeError("""Your sweep selection for kmer size
                    contained the even number or float: {0}. Please remove all even
                    numbers or floats from your input for mer_size.""".format(errors))


        #if they're supposed to be ints, convert them to ints first
        if self.sweep in self.sweep_support_int:
            self.sList = [int(x) for x in sList]

        # -------------------------Name Parameters------------------------------
          #used if-else to make this more robust
        if asPrefix is None:
            self.as_a = "as"
        else:
            self.as_a = asPrefix
        if asSI is None:
            self.as_i = 0
        else:
            self.as_i = int(asSI)
        #http://stackoverflow.com/questions/311627/
        self.as_d = tfmt("%Y%m%d")
        self.as_z = "ME"
        self.find_illegal_characters(genus, species)
        self.as_g = genus if genus else None
        self.as_s = species if species else None

        #-----------------------Directory Parameters----------------------------
        self.cwd = os.getcwd() 

    def find_illegal_characters(self, genus, species):
        #only allow the user to input ASCII letters and digits (no punctuation)
        # in the genus and species names to avoid filename errors
        illegal = [x for x in
                   str(genus).replace("None","") + str(species).replace("None","")
                   if x not in ascii_letters + digits]
        if illegal:
            raise IOError("""<genus> or <species/sample> contains the
                illegal characters: {}""".format(illegal))

    def name_gen(self, as_number, k, diploid_mode):
        """See the sweeper_output() method doctstring for format details."""
        #how to add prefixes to numbers
        #http://stackoverflow.com/a/135157
        #         ai_d_z_g_s_k_p.config

        parts =  {"1ai": "{0}{1:03d}_".format(self.as_a, as_number),
                  "2d" : "{0}_".format(self.as_d),
                  "3z" : "ME_",
                  "4g" : "{}_".format(self.as_g) if self.as_g else "",
                  "5s" : "{}_".format(self.as_s) if self.as_s else "",
                  "6k" : "k{}_".format(k),
                  "7p" : "d{}".format(diploid_mode)}
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
        p = ploid-type (diploid_mode)

        ai_d_z_g_s_k_p

        For example, if I enter the following parameters for a k=21 run on
        20160811:
          - g = "Malacosteus"
          - s = "niger"

        The assembly directory will be saved in the current working
        directory as:
          - configs/as000_20160811_ME_mala_nige_k21_p0

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
        if self.triplet:
            self.sList = [k for sublist in [[x] * 3 for x in self.sList] for k in sublist]

        self.subParams = []
        diploid_counter = 0
        asNumCounter = self.as_i
        for kmer in self.sList:
            subPDict = {}
            #set mer_size
            if self.sweep == "mer_size":
                subPDict["mer_size"] = kmer
            else:
                subPDict["mer_size"] = self.params.get("mer_size")
            #set diploid_mode
            if self.triplet:
                subPDict["diploid_mode"] = diploid_counter
            else:
                subPDict["diploid_mode"] = self.params.get("diploid_mode")
            #set assembly number
            subPDict["assem_num"] = asNumCounter
            #set assembly name
            subPDict["assem_name"] = self.name_gen(subPDict.get("assem_num"),
                                                   subPDict.get("mer_size"),
                                                   subPDict.get("diploid_mode"))
            self.subParams.append(subPDict)
            #keep track of the diploid modes to use in case of self.triplet
            if diploid_counter == 3:
                diploid_counter = 0
            else:
                diploid_counter += 1
            #increment the assembly number counter
            asNumCounter += 1

        #This block actually saves the config files
        for assemParam in self.subParams:
            with open(os.path.join(config_dir, "{0}.config".format(
                    assemParam.get("assem_name"))),
                    "w") as f:
                for key in self.params:
                    if key == "lib_seq":
                        for each in self.params["lib_seq"]:
                            #print(str(each))
                            print("lib_seq {0}".format(str(each)), file = f)
                    else:
                        if key in assemParam:
                            print("{0} {1}".format(key, assemParam.get(key)),
                                  file=f)
                        else:
                            print("{0} {1}".format(key, self.params[key]),
                                  file=f)

        # 3 return dict of {"<run string>": "<config abs path>"}.
        # THIS IS THE FINAL FUNCTIONAL OUTPUT OF THIS CLASS.
        return {name:os.path.join(config_dir, "{0}.config".format(name))
                for name in [subDict["assem_name"] for subDict in self.subParams]}

class ConfigParse:
    """This class is a Meraculous config file parser. It returns a params
    dictionary that contains all the information for one Meraculous run"""

    def __init__(self, inputFile):
        self.inputFile = inputFile
        if self.inputFile.endswith(".yaml"):
            with open(self.inputFile,'r') as infile:
                self.params = yaml.load(infile)
        else:
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
            self.config_specified = [x for x in
                                     [u for u in self.params
                                      if isinstance(self.params[u], Number)]
                                     if self.params[x] < 0]
            self.config_unspecified = [x for x in
                                       [u for u in self.params
                                        if isinstance(self.params[u], Number)]
                                       if self.params[x] >= 0] + [y for y in
                                       [v for v in self.diploid_mode
                                        if isinstance(self.diploid_mode[v], Number)]
                                       if self.diploid_mode[y] >= 0]
            self.read_config()

    def space_error(self, line):
        raise ValueError(
            """ERROR: There are too many spaces in this line in your config
            file: {0}""".format(line))

    def save_yaml(self, outFile):
        """saves the config parameters to a json file"""
        with open(outFile,'w') as myfile:
            print(yaml.dump(self.params), file=myfile)

    def sym_reads_new_config(self, newDir, sym=False, mv=False):
        """This moves the read files and renames the glob, outputs a new
        ConfigParse object with updated values"""
        newParams = copy.copy(self)
        for i in range(0,len(self.params["lib_seq"])):
            lib_seq = self.params["lib_seq"][i]
            #update the glob
            newParams.params["lib_seq"][i]["globs"] = [os.path.join(
                newDir,os.path.basename(x)) for x in lib_seq["globs"]]
            #now update the file paths and move/symlink if necessary
            for j in range(0, len(lib_seq["pairs"])):
                forward, reverse = lib_seq["pairs"][j]
                new_forward = os.path.join(newDir, os.path.basename(forward))
                new_reverse = os.path.join(newDir, os.path.basename(reverse))
                #move or symlink the files
                if sym:
                    for each in [new_forward, new_reverse]:
                        if os.path.exists(each):
                            raise ValueError("""Attempted to make a symlink for
                            a reads file, but it already exists. Chances are
                            that you used the same reads for two lines in the
                            Meraculous config file""")
                    os.symlink(forward, new_forward)
                    os.symlink(reverse, new_reverse)
                #update the reads in the new config
                newParams.params["lib_seq"][i]["pairs"][j] = [new_forward, new_reverse]
        return newParams

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
                                LS = LibSeq(line, self.inputFile)
                                #right here you need to pass the actual LS object
                                # to take adavantage of its __repr__ method
                                self.params["lib_seq"].append(LS)
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
