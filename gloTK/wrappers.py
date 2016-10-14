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

"""title: commandline_tools
authr: darrin schultz

This module:
  - defines classes that interact with command line tools

Current classes
"""
import gzip
import os
import subprocess
import sys
import time

#from itertools import chain

from gloTK import utils


#started 4:15PM
#ended at 5:30PM
#started at 6:25PM
# ended at 7:05PM
# started at 9PM
# ended at 12PM

class BaseWrapper:
    """
    This class taken from biolite (https://bitbucket.org/caseywdunn/biolite/overview)
    under GPLv3.

    A base class that handles generic wrapper functionality.

    Wrappers for specific programs should inherit this class, call `self.init`
    to specify their `name` (which is a key into the executable entries in the
    BioLite configuration file), and append their arguments to the `self.args`
    list.

    By convention, a wrapper should call `self.run()` as the final line in its
    `__init__` function. This allows for clean syntax and use of the wrapper
    directly, without assigning it to a variable name, e.g.

    wrappers.MyWrapper(arg1, arg2, ...)

    When your wrapper runs, BaseWrapper will do the following:

    * log the complete command line to diagnostics;
    * optionally call the program with a version flag (invoked with `version`)
      to obtain a version string, then log this to the :ref:`programs-table`
      along with a hash of the binary executable file;
    * append the command's stderr to a file called `name`.log in the CWD;
    * also append the command's stdout to the same log file, unless you set
      `self.stdout`, in which case stdout is redirected to a file of that name;
    * on Linux, add a memory profiling library to the LD_PRELOAD environment
      variable;
    * call the command and check its return code (which should be 0 on success,
      unless you specify a different code with `self.return_ok`), optionally
      using the CWD specified in `self.cwd` or the environment specified in
      `self.env`.
    * parse the stderr of the command to find [biolite.profile] markers and
      use the rusage values from `utils.safe_call` to populate a profile
      entity in the diagnostics with walltime, usertime, systime, mem, and
      vmem attributes.
    """

    def __init__(self, name, **kwargs):
        self.name = name
        self.shell = '/bin/sh'
        self.args = []
        self.return_ok = kwargs.get('return_ok', 0)
        self.cwd = kwargs.get('cwd', os.getcwd())
        self.stdout = kwargs.get('stdout')
        self.outdir = kwargs.get('outDir')
        self.gzip = kwargs.get('gzip', None)
        print("init")
        self.stdout_append = kwargs.get('stdout_append')
        self.pipe = kwargs.get('pipe')
        self.env = os.environ.copy()
        self.max_concurrency = kwargs.get('max_concurrency', 1)


    init = __init__
    """A shortcut for calling the BaseWrapper __init__ from a subclass."""

    def check_input(self, flag, path):
        """
        turns path into an absolute path and checks that it exists, then
        appends it to the args using the given flag (or none).
        """
        path = os.path.abspath(path)
        if os.path.exists(path):
            if flag:
                self.args.append(flag)
            self.args.append(path)
        else:
            utils.die("input file for flag '%s' does not exists:\n  %s" % (
                      flag, path))

    def check_path(self, path):
        """
        turns path into an absolute path and checks that it exists, then
        returns it as a string.
        """
        path = os.path.abspath(path)
        if os.path.exists(path):
            return path
        else:
            utils.die("input file does not exists:\n  {}".format(path))


    # def add_threading(self, flag):
    #     """
    #     Indicates that this wrapper should use threading by appending an
    #     argument with the specified `flag` followed by the number of threads
    #     specified in the BioLite configuration file.
    #     """
    #     threads = min(int(config.get_resource('threads')), self.max_concurrency)
    #     if threads > 1:
    #         self.args.append(flag)
    #         self.args.append(threads)


    # def add_openmp(self):
    #     """
    #     Indicates that this wrapper should use OpenMP by setting the
    #     $OMP_NUM_THREADS environment variable equal to the number of threads
    #     specified in the BioLite configuration file.
    #     """
    #     threads = min(int(config.get_resource('threads')), self.max_concurrency)
    #     self.env['OMP_NUM_THREADS'] = str(threads)

    def run(self):
        """
        Call this function at the end of your class's `__init__` function.
        """
        stderr = os.path.abspath(os.path.join(self.outdir, self.name + '.log'))

        if self.pipe:
            self.args += ('|', self.pipe, '2>>'+stderr)

        if self.gzip:
            self.args += ('|', 'gzip', '1>', self.gzip)
        else:
            self.args.append('2>>'+stderr)
            self.args.append('1>>'+stderr)

        # Print timestamp to log

        log = open(stderr, 'a')
        log.write("[gloTK] timestamp={}\n".format(utils.timestamp()))

        cmd = ' '.join(map(str, self.args))
        print(cmd)
        log.write(cmd)

        start = time.time()
        save_cwd = os.getcwd()
        try:
            utils.safe_mkdir(self.outdir)
            os.chdir(self.outdir)
            spawn_pid = os.spawnle(os.P_NOWAIT, self.shell, self.shell, '-c', cmd, self.env)
            wait_pid, retcode, rusage = os.wait4(spawn_pid, 0)
            if wait_pid != spawn_pid:
                utils.die("could not wait for process %d: got %d" % (spawn_pid, wait_pid))
            os.chdir(save_cwd)
        except OSError as e:
            utils.info(e)
            utils.die("could not run wrapper for command:\n%s" % cmd)

        elapsed = time.time() - start
        retcode = os.WEXITSTATUS(retcode)

        if (self.return_ok is not None) and (self.return_ok != retcode):
            # Give some context to the non-zero return.
            if os.path.isfile(stderr):
                subprocess.call(['tail', '-3', stderr])
            utils.die("non-zero return (%d) from command:\n%s" % (retcode, cmd))
        log.close()

    def version_jar(self):
        """
        Special case of version() when the executable is a JAR file.
        """
        cmd = config.get_command('java')
        cmd.append('-jar')
        #self.version(cmd=cmd, path=self.cmd[0])


    def run_jar(self, mem=None):
        """
        Special case of run() when the executable is a JAR file.
        """
        cmd = config.get_command('java')
        if mem:
            cmd.append('-Xmx%s' % mem)
        cmd.append('-jar')
        cmd += self.cmd
        self.run(cmd)


class Seqprep (BaseWrapper):
    """
    usage: Seqprep(required_args, other_args)

    This doctring is modified from SeqPrep2 unix command line -h option.

    Required Arguments:
      forwardPath    <first read input fastq filepath>
      reversePath    <second read input fastq filepath>
      forwardOutFile <first read output fastq filename>
      reverseOutFile <second read output fastq filename>
      outDir         <directory where files will be saved

    Arguments for Adapter/Primer Trimming (Optional):
      qualCutoff     <Quality score cutoff for mismatches to be counted in
                       overlap; default = 13>
      lenCutoff      <Minimum length of a trimmed or merged read to print it;
                       default = 30>
      forwardAdapter <forward read primer/adapter sequence to trim as it would
                       appear at the end of a read (recommend about 20bp of this)
                       (should validate by grepping a file); default (genomic
                       non-multiplexed adapter1) = "AGATCGGAAGAGCACACGTC">
      reverseAdapter <reverse read primer/adapter sequence to trim as it would
                       appear at the end of a read (recommend about 20bp of this)
                       (should validate by grepping a file); default (genomic
                       non-multiplexed adapter2) = "AGATCGGAAGAGCGTCGTGT">

    Optional Arguments for Merging
      mergedOutFile  <perform merging and output merged reads to this file>
      overlapMin     <minimum overall base pair overlap to merge two reads;
                       default = 15>
      prettyOutFile  <write pretty alignments to this file for visual
                       Examination>


    Optional Arguments for Rejecting Unmerged Sequences with Adapter Primers
      editReject     <perform sequence match rejection and set cutoff edit distance
                       for rejection sequences; default = 1>
      forwardReject  <first read primer rejection sequence;
                       default = "AGATCGGAAGAGCACACGTC">
      reverseReject  <second read primer rejection sequence;
                       default = "AGATCGGAAGAGCGTCGTGT">
    """
    def __init__(self, **kwargs):
        self.init('seqprep2', **kwargs)
        for each in kwargs:
            if "Path" in each:
                kwargs[each] = self.check_path(kwargs[each])
            if "File" in each:
                kwargs[each] = os.path.join(kwargs["outDir"], kwargs[each])

        self.args = ["SeqPrep2",
                "-f", kwargs["forwardPath"],
                "-r", kwargs["reversePath"],
                "-1", kwargs["forwardOutFile"],
                "-2", kwargs["reverseOutFile"],
                "-q", kwargs.get("qualCutoff", 13),
                "-L", kwargs.get("lenCutoff", 30),
                "-A", kwargs.get("forAdapter", "AGATCGGAAGAGCACACGTC"),
                "-B", kwargs.get("revAdapter", "AGATCGGAAGAGCGTCGTGT"),
                "-d", kwargs.get("editReject", 1),
                "-C", kwargs.get("forwardReject", "AGATCGGAAGAGCACACGTC"),
                "-D", kwargs.get("reverseReject", "AGATCGGAAGAGCGTCGTGT")]

        if kwargs.get("mergedOutFile"):
            self.args += ["-s", kwargs["mergedOutFile"],
                          "-E", kwargs["prettyOutFile"],
                          "-x", kwargs.get("prettyNum", 50),
                          "-o", kwargs.get("overlapMin", 30)]

        self.run()


class Seqtk(BaseWrapper):
    """
    Required Arguments:
      inputPath
      outDir
      readCount

    Optional Arguments:
      seed
    """
    def __init__(self, **kwargs):
        """Run seqtk to sample the reads and return output if there is any"""
        self.init('seqtk', **kwargs)

        print("instance")
        self.gzip = os.path.join(kwargs["outDir"],
            "{}_{}reads.fastq.gz".format(
                utils.fastx_basename(kwargs["inputPath"]),
                kwargs["readCount"]) )
        print(self.gzip)

        try:
            os.remove(self.gzip)
        except OSError:
            pass

        self.args = ['seqtk', 'sample',
                     '-s', kwargs.get("seed", 100),
                     kwargs["inputPath"],
                     kwargs["readCount"]]
        self.run()
