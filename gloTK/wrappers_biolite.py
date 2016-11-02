# BioLite - Tools for processing gene sequence data and automating workflows
# Copyright (c) 2012-2014 Brown University. All rights reserved.
#
# This file is part of BioLite.
#
# BioLite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BioLite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioLite.  If not, see <http://www.gnu.org/licenses/>.

"""
A series of wrappers for external calls to various bioinformatics tools.
"""

import glob
import math
import operator
import os
import random
import shlex
import subprocess
import sys
import time

from collections import namedtuple
from itertools import chain

import config
import diagnostics
import utils

class BaseWrapper:
	"""
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
		self.cmd = config.get_command(name)
		self.args = []
		self.return_ok = kwargs.get('return_ok', 0)
		self.cwd = kwargs.get('cwd', os.getcwd())
		self.stdout = kwargs.get('stdout')
		self.stdout_append = kwargs.get('stdout_append')
		self.pipe = kwargs.get('pipe')
		self.env = os.environ.copy()
		self.max_concurrency = kwargs.get('max_concurrency', sys.maxint)
		self.output_patterns = None


	init = __init__
	"""A shortcut for calling the BaseWrapper __init__ from a subclass."""


	def check_input(self, flag, path):
		"""
		Turns path into an absolute path and checks that it exists, then
		appends it to the args using the given flag (or None).
		"""
		path = os.path.abspath(path)
		if os.path.exists(path):
			if flag:
				self.args.append(flag)
			self.args.append(path)
		else:
			utils.die("input file for flag '%s' does not exists:\n  %s" % (
			          flag, path))


	def add_threading(self, flag):
		"""
		Indicates that this wrapper should use threading by appending an
		argument with the specified `flag` followed by the number of threads
		specified in the BioLite configuration file.
		"""
		threads = min(int(config.get_resource('threads')), self.max_concurrency)
		if threads > 1:
			self.args.append(flag)
			self.args.append(threads)


	def add_openmp(self):
		"""
		Indicates that this wrapper should use OpenMP by setting the
		$OMP_NUM_THREADS environment variable equal to the number of threads
		specified in the BioLite configuration file.
		"""
		threads = min(int(config.get_resource('threads')), self.max_concurrency)
		self.env['OMP_NUM_THREADS'] = str(threads)


	def version(self, flag=None, cmd=None, path=None):
		"""
		Generates and logs a hash to distinguish this particular installation
		of the program (on a certain host, with a certain compiler, program
		version, etc.)

		Specify the optional 'binary' argument if the wrapper name is not
		actually the program, e.g. if your program has a Perl wrapper script.
		Set 'binary' to the binary program that is likely to change between
		versions.

		Specify the optional 'cmd' argument if the command to run for version
		information is different than what will be invoked by `run` (e.g.
		if the program has a perl wrapper script, but you want to version an
		underlying binary executable).
		"""

		# Setup the command to run.
		if not cmd:
			cmd = list(self.cmd)
		if flag:
			cmd.append(flag)

		# Run the command.
		try:
			vstring = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
		except subprocess.CalledProcessError as e:
			vstring = e.output
		except OSError as e:
			utils.failed_executable(cmd[0], e)

		if not path:
			path = cmd[0]

		# Generate a hash.
		vhash = diagnostics.log_program_version(self.name, vstring, path)
		if vhash:
			diagnostics.prefix.append(self.name)
			diagnostics.log('version', vhash)
			diagnostics.prefix.pop()


	def version_jar(self):
		"""
		Special case of version() when the executable is a JAR file.
		"""
		cmd = config.get_command('java')
		cmd.append('-jar')
		cmd += self.cmd
		self.version(cmd=cmd, path=self.cmd[0])


	def run(self, cmd=None):
		"""
		Call this function at the end of your class's `__init__` function.
		"""
		diagnostics.prefix.append(self.name)

		if not cmd:
			cmd = self.cmd

		stderr = os.path.abspath(self.name + '.log')
		self.args.append('2>>'+stderr)

		if self.pipe:
			self.args += ('|', self.pipe, '2>>'+stderr)

		# Write to a stdout file if it was set by the derived class.
		# Otherwise, stdout and stderr will be combined into the log file.
		if self.stdout:
			stdout = os.path.abspath(self.stdout)
			self.args.append('1>'+stdout)
			diagnostics.log('stdout', stdout)
		elif self.stdout_append:
			stdout = os.path.abspath(self.stdout_append)
			self.args.append('1>>'+stdout)
			diagnostics.log('stdout', stdout)
		else:
			self.args.append('1>>'+stderr)

		# Print timestamp to log
		open(stderr, 'a').write("[biolite] timestamp=%s\n" % utils.timestamp())
		diagnostics.log('log', stderr)

		cmd = ' '.join(chain(cmd, map(str, self.args)))
		diagnostics.log('command', cmd)

		start = time.time()
		save_cwd = os.getcwd()
		try:
			os.chdir(self.cwd)
			spawn_pid = os.spawnle(os.P_NOWAIT, self.shell, self.shell, '-c', cmd, self.env)
			wait_pid, retcode, rusage = os.wait4(spawn_pid, 0)
			if wait_pid != spawn_pid:
				utils.die("could not wait for process %d: got %d" % (spawn_pid, wait_pid))
			os.chdir(save_cwd)
		except OSError as e:
			utils.info(e)
			utils.die("could not run wrapper for command:\n%s" % cmd)
			#utils.failed_executable(exe, e)

		elapsed = time.time() - start
		retcode = os.WEXITSTATUS(retcode)

		if (self.return_ok is not None) and (self.return_ok != retcode):
			# Give some context to the non-zero return.
			if os.path.isfile(stderr):
				subprocess.call(['tail', '-3', stderr])
			utils.die("non-zero return (%d) from command:\n%s" % (retcode, cmd))

		# Log profile.
		diagnostics.prefix.append('profile')
		diagnostics.log('name', self.name)
		diagnostics.log('return', retcode)
		diagnostics.log('walltime', elapsed)
		diagnostics.log('usertime', rusage.ru_utime)
		diagnostics.log('systime', rusage.ru_stime)
		if config.uname == 'Darwin':
			diagnostics.log('maxrss', rusage.ru_maxrss / 1024)
		else:
			diagnostics.log('maxrss', rusage.ru_maxrss)
		diagnostics.prefix.pop()

		# Reverse any output patterns, since they will be matched against
		# program output from the last line backward.
		if self.output_patterns:
			self.output_patterns.reverse()

		diagnostics.log_program_output(stderr, self.output_patterns)
		diagnostics.prefix.pop()


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


### BioLite command line tools ###


class Coverage (BaseWrapper):
	"""
	usage: coverage [-i SAM] [-o STATS]

	Parses a SAM alignment file and writes a coverage table to STATS with
	columns for the reference name, the length of the referene, and the number
	of reads covering it in the alignment.
	"""

	def __init__(self, input, *args, **kwargs):
		self.init('coverage', **kwargs)
		self.version('-v')
		self.check_input('-i', input)
		self.args += args
		self.run()


class Exclude (BaseWrapper):
	"""
	usage: exclude -x EXCLUDE_FILE [-k] [...] [-i INPUT ...] [-o OUTPUT ...]

	Filters all the reads in the input files (FASTA or FASTQ is automatically
	detected) and excludes those with ids found in any of the EXCLUDE_FILEs.

	If multiple input files are specified, these are treated as paired files.
	So if a sequence in one input is excluded, its pair is also excluded from
	the same position in all other input files.

	If the -k flag is specified, invert the selection to keep instead of exclude.
	"""

	def __init__(self, excludes, inputs, outputs, *args, **kwargs):
		self.init('exclude', **kwargs)
		self.version('-v')
		for x in excludes:
			self.check_input('-x', x)
		for i in inputs:
			self.check_input('-i', i)
		for o in outputs:
			self.args += ('-o', o)
		self.args += args
		self.run()


class Fastq2Fasta (BaseWrapper):
	"""
	usage: fastq2fasta -i FASTQ [...] [-o FASTA ...] [-q QUAL ...] [-a]
	                   [-t OFFSET] [-s SUFFIX]

	Converts each FASTQ input file to a FASTA file and quality score file
	with the names <basename>.fasta and <basename>.fasta.qual, where <basename>
	is the name of INPUT up to the last period (or with the names FASTA and QUAL
	if specified).

	FASTA and QUAL are *appended* to (not truncated).
	"""

	def __init__(self, input, *args, **kwargs):
		self.init('fastq2fasta', **kwargs)
		self.version('-v')
		self.check_input('-i', input)
		self.args += args
		self.run()


class Fasta2Fastq (BaseWrapper):
	"""
	usage: fasta2fastq -i FASTA [...] -q QUAL [...] [-o FASTQ] [-a] [-t OFFSET]

	Merges each FASTA file and its corresponding QUAL file into a FASTQ file
	with the name <basename>.fastq, where <basename> in the FASTA name up to the
	last period (or with name FASTQ if specified).

	FASTQ is *appended* to (not truncated).
	"""

	def __init__(self, input, qual, *args, **kwargs):
		self.init('fasta2fastq', **kwargs)
		self.version('-v')
		self.check_input('-i', input)
		self.check_input('-q', qual)
		self.args += args
		self.run()


class FilterIllumina (BaseWrapper):
	"""
	usage: filter_illumina [-i INPUT ...] [-o OUTPUT ...] [-u UNPAIRED-OUTPUT] [-f]
	                       [-t OFFSET] [-q QUALITY] [-n NREADS] [-a] [-b] [-s SEP]

	Filters out low-quality and adapter-contaminated reads from Illumina data.

	If multiple input files are specified, these are treated as paired files.
	So if a sequence in one input is filtered, its pair is also filtered from
	the same position in all other input files.
	"""

	def __init__(self, inputs, outputs, *args, **kwargs):
		self.init('filter_illumina', **kwargs)
		self.version('-v')
		for i in inputs:
			self.check_input('-i', i)
		for o in outputs:
			self.args += ('-o', o)
		self.args += args
		self.run()


class InsertStats (BaseWrapper):
	"""
	usage: insert_stats -i SAM -o HIST -m MAX_INSERT

	Reads a SAM alignment file and uses it to estimate the mean and std. dev.
	of the insert size of the mapped paired-end reads. A histogram of all insert
	sizes encountered is written to the HIST file.
	"""

	def __init__(self, input, *args, **kwargs):
		self.init('insert_stats', **kwargs)
		self.version('-v')
		self.check_input('-i', input)
		self.args += ('-m', kwargs.get('max_insert',
		                               config.get_resource('max_insert_size')))
		self.args += args
		self.run()


class Interleave (BaseWrapper):
	"""
	usage: interleave -i INPUT [...] [-o OUTPUT] [-s SEP]

	Interleaves the records in the input files (FASTA or FASTQ is automatically
	detected) and writes them to OUTPUT, or to stdout if no OUTPUT is specified.
	"""

	def __init__(self, inputs, output, *args, **kwargs):
		self.init('interleave', **kwargs)
		self.version('-v')
		for i in inputs:
			self.args += ('-i', i)
		self.args += ('-o', output)
		self.args += args
		self.run()


class Randomize (BaseWrapper):
	"""
	usage: randomize [-i INPUT] [-o OUTPUT] [-r READ-ORDER] [-w WRITE-ORDER]

	Randomizes the order of sequences in each INPUT file and writes these to a
	corresponding OUTPUT file. By default, a new random write order is generated
	and saved to WRITE-ORDER, if specified. Alternatively, specifying a
	READ-ORDER file uses that order instead of a random one.
	"""

	def __init__(self, input, *args, **kwargs):
		self.init('randomize', **kwargs)
		self.version('-v')
		self.check_input('-i', input)
		self.args += args
		self.run()


### Third-party command line tools ###


class Abacas (BaseWrapper):
	"""
	ABACAS: Algorithm Based Automatic Contiguation of Assembled Sequences
	http://abacas.sourceforge.net
	"""

	def __init__(self, contigs, reference, program, *args, **kwargs):
		self.init('abacas', **kwargs)
		self.args += ('-q', contigs, '-r', reference, '-p', program) + args
		self.run()


class BlastN (BaseWrapper):
	"""
	blastn from NCBI Blast+
	http://blast.ncbi.nlm.nih.gov/
	"""

	def __init__(self, query, db, *args, **kwargs):
		self.init('blastn', **kwargs)
		self.version('-version')
		self.args += ('-query', query, '-db', os.path.abspath(db))
		self.add_threading('-num_threads')
		self.args += args
		self.run()


class BlastP (BaseWrapper):
	"""
	blastp from NCBI Blast+
	http://blast.ncbi.nlm.nih.gov/
	"""
	def __init__(self, query, db, *args, **kwargs):
		self.init('blastp', **kwargs)
		self.version('-version')
		self.args += ('-query', query, '-db', os.path.abspath(db))
		self.add_threading('-num_threads')
		self.args += args
		self.run()


class BlastX (BaseWrapper):
	"""
	blastx from NCBI Blast+
	http://blast.ncbi.nlm.nih.gov/
	"""

	def __init__(self, query, db, *args, **kwargs):
		self.init('blastx', **kwargs)
		self.version('-version')
		self.args += ('-query', query, '-db', os.path.abspath(db))
		self.add_threading('-num_threads')
		self.args += args
		self.run()


class Bowtie2 (BaseWrapper):
	"""
	A wrapper for the bowtie2 short-read aligner.
	http://bowtie-bio.sourceforge.net/

	For paired inputs, you can specify the maximum insert size (e.g. the
	length of the gap between the reads) with the 'max_insert' keyword
	argument. If you don't specify one, the diagnostics database will be
	searched for a previous run of the 'insert_size' pipeline for an estimate.
	"""

	def __init__(self, inputs, db, *args, **kwargs):
		self.init('bowtie2', **kwargs)
		self.version('--version', config.get_command('bowtie2-align'))
		if isinstance(inputs, basestring):
			self.check_input('-U', inputs)
		elif len(inputs) == 1:
			self.check_input('-U', inputs[0])
		elif len(inputs) == 2:
			self.check_input('-1', inputs[0])
			self.check_input('-2', inputs[1])
			self.args.append('-X')
			self.args.append(
				kwargs.get('max_insert', diagnostics.lookup_insert_size().max))
		else:
			utils.die("Bowtie2 wrapper expects either 1 (SE) or 2 (PE) inputs")
		self.args += ('-x', db)
		self.add_threading('-p')
		self.args += args
		self.output_patterns = map(diagnostics.OutputPattern._make, [
			(r"(\d+) reads; of these:$",0,"nreads"),
			(r"  (\d+) \S+ were paired; of these:$",0,"npairs"),
			(r"    (\d+) \S+ aligned concordantly 0 times$",0,"nconcord0"),
			(r"    (\d+) \S+ aligned concordantly exactly 1 time$",0,"nconcord1"),
			(r"    (\d+) \S+ aligned concordantly >1 times$",0,"nconcord2"),
			(r"      (\d+) \S+ aligned discordantly 1 time$",0,"ndiscord1"),
			(r"      (\d+) mates make up the pairs; of these:$",0,"nunpaired"),
			(r"        (\d+) \S+ aligned 0 times$",0,"nunpaired0"),
			(r"        (\d+) \S+ aligned exactly 1 time$",0,"nunpaired1"),
			(r"        (\d+) \S+ aligned >1 times$",0,"nunpaired2")])
		self.run()


class Bowtie2Build (BaseWrapper):
	"""
	A wrapper for bowtie2-build component of Bowtie2.
	http://bowtie-bio.sourceforge.net/
	"""

	def __init__(self, input, db, *args, **kwargs):
		self.init('bowtie2-build', **kwargs)
		self.version('--version')
		self.check_input(None, input)
		self.args.append(db)
		self.args += args
		self.run()


class Chrysalis (BaseWrapper):
	"""
	The Chrysalis component of the Trinity RNA-seq assembler:
	http://trinityrnaseq.sourceforge.net
	"""

	def __init__(self, input, iworm, *args, **kwargs):
		self.init('chrysalis', **kwargs)
		self.add_threading('-cpu')
		self.add_openmp()
		self.check_input('-i', input)
		self.check_input('-iworm', iworm)
		mem = utils.mem_to_mb(config.get_resource('memory'))
		self.args += ('-sort_buffer_size', '%dM' % int(0.8 * mem))
		self.args += args
		self.run()


class Dustmasker (BaseWrapper):
	"""
	dustmasker from NCBI Blast+
	http://blast.ncbi.nlm.nih.gov/
	"""

	def __init__(self, input, *args, **kwargs):
		self.init('dustmasker', **kwargs)
		self.version('-version-full')
		self.check_input('-in', input)
		self.args += args
		self.run()


class FastQC (BaseWrapper):
	"""
	A wrapper for FastQC.
	http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/
	"""

	def __init__(self, input, *args, **kwargs):
		self.init('fastqc', **kwargs)
		self.version('-v')
		self.add_threading('-t')
		self.args += args
		self.check_input(None, input)
		self.run()


class Gblocks (BaseWrapper):
	"""
	Selection of conserved block from multiple sequence alignments for
	phylogenetics
	http://molevol.cmima.csic.es/castresana/Gblocks.html
	"""

	def __init__(self, input, *args, **kwargs):
		self.init('Gblocks', **kwargs)
		# Ignore Gblocks broken exit code
		self.return_ok = None
		self.check_input(None, input)
		self.args += args
		self.run()


class Inchworm (BaseWrapper):
	"""
	The inchworm component of the Trinity RNA-seq assembler:
	http://trinityrnaseq.sourceforge.net
	"""

	def __init__(self, mode, input, *args, **kwargs):
		self.init('inchworm', **kwargs)
		self.check_input(mode, input)
		self.add_openmp()
		self.args.append('--run_inchworm')
		self.args += args
		self.run()


class JellyfishCount (BaseWrapper):
	"""
	"""

	def __init__(self, input, kmer, *args, **kwargs):
		self.init('jellyfish', **kwargs)
		self.version('--version')
		#
		# From the Jellyfish manual (section 2.3), hash size in bytes is:
		#
		#   2^l * (2k-l+r+1)/8
		#
		# s = 2^l is the hash size parameter given by -s. By default, r=5 and
		# letting l be 0:
		#
		#   s ~= 8 * mem / (2k + 6)
		#
		mem = 1048576 * utils.mem_to_mb(config.get_resource('memory'))
		mem = 2 ** int(math.log(8 * mem / (2*kmer + 6) , 2))
		self.args += ('count', '-m', kmer, '-s', mem)
		self.add_threading('-t')
		self.check_input(None, input)
		self.args += args
		self.run()


class JellyfishDump (BaseWrapper):
	"""
	"""

	def __init__(self, input, *args, **kwargs):
		self.init('jellyfish', **kwargs)
		self.version('--version')
		self.args.append('dump')
		self.check_input(None, input)
		self.args += args
		self.run()


class MultiBlast (BaseWrapper):
	"""
	usage: multiblast BLAST THREADS QUERY_LIST OUT [ARGS]

	Runs a Blast PROGRAM (e.g. blastx, blastn, blastp) in parallel on a list of
	queries (in QUERY_LIST).  Additional arguments to PROGRAM can be appended as
	ARGS.

	The PROGRAM is called on each query with threading equal to THREADS.
	Recommendation: set THREADS to the number of cores divided by the number of
	query files.

	The individual XML outputs for each query file are concatenated into a single
	output file OUT.

	Example usage:
	multiblast blastn 4 "query1.fa query2.fa" all-queries.xml -e 1e-6
	"""

	def __init__(self, blast, threads, qlist, db, out, evalue=0.0001, targets=20):
		if not glob.glob(db + '.*'):
			utils.die("missing blast database '%s'" % db)
		self.init(blast)
		self.version('-version')
		self.args += (
			threads, ' '.join(qlist), out,
			'-db', os.path.abspath(db), '-evalue', evalue,
			'-max_target_seqs', targets)
		self.run(config.get_command('multiblast') + self.cmd)


class Macse (BaseWrapper):
	"""
	Multiple alignment of coding sequences.
	"""

	def __init__(self, input, output, *args, **kwargs):
		self.init('macse')
		self.version_jar()
		self.check_input('-i', input)
		self.args += ('-o', output)
		self.args += args
		mem = 0.9 * utils.mem_to_mb(config.get_resource('memory'))
		self.run_jar('%dM' % mem)


class MakeBlastDB (BaseWrapper):
	"""
	makeblastdb from NCBI Blast+
	http://blast.ncbi.nlm.nih.gov/
	"""

	def __init__(self, input, db, dbtype, *args, **kwargs):
		self.init('makeblastdb', **kwargs)
		self.version('-version')
		self.check_input('-in', input)
		self.args += ('-dbtype', dbtype, '-out', os.path.abspath(db))
		self.args += args
		self.run()


class Mcl (BaseWrapper):
	"""
	Markov Clustering Algorithm (MCL) for analysis of networks
	http://micans.org/mcl/
	"""

	def __init__(self, input, *args, **kwargs):
		self.init('mcl', **kwargs)
		self.version('--version')
		self.check_input(None, input)
		self.add_threading('-te')
		self.args += args
		self.run()


class Minimo (BaseWrapper):
	"""
	Minimo: overlap graph assembler for small data sets
	From the AMOS assembler package.
	http://amos.sourceforge.net
	"""

	def __init__(self, fasta, *args, **kwargs):
		self.init('Minimo', **kwargs)
		self.args.append(fasta)
		self.args += args
		self.run()


class Oases (BaseWrapper):
	"""
	Oases, a *de novo* transcriptome assembler
	http://www.ebi.ac.uk/~zerbino/oases/
	"""

	def __init__(self, workdir, *args, **kwargs):
		self.init('oases')
		self.version()
		self.args.append(workdir)
		self.args += args
		self.run()


class Oma (BaseWrapper):
	"""
	"""

	def __init__(self, **kwargs):
		self.init('oma', **kwargs)
		parameters = """
			ReuseCachedResults := true;
			NP := %d;
			MinScore := 181;
			LengthTol := 0.61;
			StablePairTol := 1.81;
			VerifiedPairTol := 1.53;
			MinSeqLen := 50;
			StableIdsForGroups := false;
			DoHierarchicalGroups := true;
			MaxTimePerLevel := 1200;
			SpeciesTree := 'estimate';
			ReachabilityCutoff := 0.65;
			UseEsprit := false;
			DistConfLevel := 2;
			MinProbContig := 0.4;
			MaxContigOverlap := 5;
			MinSeqLenContig := 20;
			MinBestScore := 250;
		""" % int(config.get_resource('threads'))
		open(os.path.join(self.cwd, 'parameters.drw'), 'w').write(parameters)
		self.run()


class PartitionChrysalis (BaseWrapper):
	"""
	The partitioning script for the Chrysalis component of the Trinity RNA-seq
	assembler:
	http://trinityrnaseq.sourceforge.net
	"""

	def __init__(self, debruijn, reads, *args, **kwargs):
		self.init('partition_chrysalis', **kwargs)
		self.check_input('--deBruijns', debruijn)
		self.check_input('--componentReads', reads)
		self.args += args
		self.run()


class Parallel (BaseWrapper):
	"""
	GNU parallel utility
	http://www.gnu.org/software/parallel/
	"""

	def __init__(self, commands, *args, **kwargs):
		self.init('parallel', **kwargs)
		self.version('--version')
		self.args += (
			'--gnu', '-a', commands, '-j',
			kwargs.get('threads', config.get_resource('threads')))
		hostlist = config.get_resource_default('hostlist', None)
		if hostlist:
			self.args += ('-S', hostlist)
			if self.cwd:
				self.args += ('--wd', self.cwd)
			else:
				self.args += ('--wd', os.getcwd())
		self.args += args
		self.run()


class Raxml (BaseWrapper):
	"""
	Maximum Likelihood based inference of phylogenetic trees.
	"""

	def __init__(self, input, *args, **kwargs):
		self.init('raxml', **kwargs)
		self.version('-v')
		self.check_input('-s', input)
		self.add_threading('-T')
		self.args += args
		self.run()


class RaxmlMpi (BaseWrapper):
	"""
	Maximum Likelihood based inference of phylogenetic trees
	(MPI version).
	"""

	def __init__(self, mpirun, input, *args, **kwargs):
		self.init('raxml-mpi', **kwargs)
		self.cmd.insert(0, mpirun)
		self.version('-v')
		self.check_input('-s', input)
		self.args += args
		self.run()


class RaxmlHybrid (BaseWrapper):
	"""
	Maximum Likelihood based inference of phylogenetic trees
	(MPI-hybrid version).
	"""

	def __init__(self, mpirun, input, *args, **kwargs):
		self.init('raxml-hybrid', **kwargs)
		self.cmd.insert(0, mpirun)
		self.version('-v')
		self.check_input('-s', input)
		self.add_threading('-T')
		self.args += args
		self.run()


class RpsBlast (BaseWrapper):
	"""
	rpsblast from NCBI Blast+
	http://blast.ncbi.nlm.nih.gov/
	"""

	def __init__(self, query, db, *args, **kwargs):
		self.init('rpsblast', **kwargs)
		self.version('-version')
		self.args += ('-query', query, '-db', os.path.abspath(db))
		self.add_threading('-num_threads')
		self.args += args
		self.run()


class RsemReference (BaseWrapper):
	"""
	http://deweylab.biostat.wisc.edu/rsem/
	"""

	def __init__(self, input, prefix, *args, **kwargs):
		self.init('rsem-prepare-reference', **kwargs)
		self.args += args
		self.args += (input, prefix)
		self.run()


class RsemExpression (BaseWrapper):
	"""
	http://deweylab.biostat.wisc.edu/rsem/
	"""

	def __init__(self, inputs, prefix, name, *args, **kwargs):
		self.init('rsem-calculate-expression', **kwargs)
		self.add_threading('--num-threads')
		max_insert = kwargs.get('max_insert')
		if len(inputs) == 2:
			if max_insert is None:
				max_insert = diagnostics.lookup_insert_size().max
			self.args += ('--paired-end', '--fragment-length-max', int(max_insert))
		self.args += args
		for input in inputs:
			if input.endswith('.gz'):
				self.shell = '/bin/bash'
				self.args.append('<(gzip -dc %s)' % input)
			else:
				self.args.append(input)
		self.args += (prefix, name)
		self.output_patterns = map(diagnostics.OutputPattern._make, [
			(r"# reads processed: (\d+)$",0,"nreads"),
			(r"# reads with at least one reported alignment: (\d+) \S+$",0,"naligned"),
			(r"# reads that failed to align: (\d+) \S+$",0,"nfailed")])
		self.run()


class SamTools (BaseWrapper):
	def __init__(self, input, *args, **kwargs):
		self.init('samtools', **kwargs)
		self.version()
		self.args += args
		self.check_input(None, input)
		self.run()


class SamView (BaseWrapper):
	def __init__(self, input_path, regions, output_path):
		self.init('samtools')
		self.version()
		self.args += ('view', '-o', output_path, input_path)
		self.args += regions
		self.run()


class SamToolsSort (BaseWrapper):
	def __init__(self, input, prefix, *args, **kwargs):
		self.init('samtools', **kwargs)
		self.version()
		self.args.append('sort')
		self.args += args
		self.check_input(None, input)
		self.args.append(prefix)
		self.run()


class SamIndex (BaseWrapper):
	def __init__(self, input_path):
		self.init('samtools')
		self.version()
		self.args += ('index', input_path)
		self.run()


class SamPileup (BaseWrapper):
	def __init__(self, reference_path, bam_path, output_path):
		self.init('samtools')
		self.version()
		self.args += (
			'mpileup', '-BQ0', '-d1000000000', '-f', reference_path, bam_path)
		self.stdout = output_path
		self.run()


class Spades (BaseWrapper):
	"""
	SPAdes de novo assembler
	http://bioinf.spbau.ru/spades
	"""
	def __init__(self, inputs, *args, **kwargs):
		self.init("spades.py", **kwargs)
		self.name = "spades"
		self.version()
		self.add_threading("-t")
		mem = 0.9 * utils.mem_to_mb(config.get_resource('memory')) / 1024
		self.args += ("-m", max(1, int(mem)))
		# Detect inputs
		if isinstance(inputs, basestring):
			self.check_input('-s', inputs)
		elif len(inputs) == 1:
			self.check_input('-s', inputs[0])
		elif len(inputs) == 2:
			self.check_input('-1', inputs[0])
			self.check_input('-2', inputs[1])
		else:
			utils.die("expected either 1 (SE) or 2 (PE) inputs")
		self.args += args
		self.run()


class Sqlite3 (BaseWrapper):
	def __init__(self, dbpath, sql, *args, **kwargs):
		self.init('sqlite3', **kwargs)
		self.version('-version')
		self.args += args
		self.check_input(None, dbpath)
		self.args.append('"%s"' % sql.replace('"', '\"'))
		self.run()


class TBlastX (BaseWrapper):
	"""
	tblastx from NCBI Blast+
	http://blast.ncbi.nlm.nih.gov/
	"""

	def __init__(self, query, db, *args, **kwargs):
		self.init('tblastx', **kwargs)
		self.version('-version')
		self.args += ('-query', query, '-db', os.path.abspath(db))
		self.add_threading('-num_threads')
		self.args += args
		self.run()


class Transdecoder (BaseWrapper):
	"""
	Identification of candidate coding sequences
	http://transdecoder.sourceforge.net
	"""

	def __init__(self, input, **kwargs):
		self.init('transdecoder', **kwargs)
		self.check_input('-t', input)
		self.run()


class Trinity (BaseWrapper):
	"""
	Trinity RNA-seq assembler
	http://trinityrnaseq.sourceforge.net
	"""

	def __init__(self, inputs, *args, **kwargs):
		self.init('trinity', **kwargs)
		self.version('--version')

		# Detect inputs
		if isinstance(inputs, basestring):
			self.check_input('--single', inputs)
		elif len(inputs) == 1:
			self.check_input('--single', inputs[0])
		elif len(inputs) == 2:
			self.check_input('--left', inputs[0])
			self.check_input('--right', inputs[1])
		else:
			utils.die("expected either 1 (SE) or 2 (PE) inputs")

		# Detect file type
		ext = os.path.splitext(self.args[-1])[1]
		if ext == '.fa':
			self.args += ('--seqType', 'fa')
		elif ext == '.fq':
			self.args += ('--seqType', 'fq')
		else:
			utils.info("warning: could not determine sequence type of inputs")

		# Java uses roughly 2 CPUs per Butterfly call with GC etc. so reduce
		# the number of threads by half.
		#threads = kwargs.get('threads', int(config.get_resource('threads')))
		#self.max_concurrency = max(1, threads/2)

		mem = utils.mem_to_mb(config.get_resource('memory'))
		self.args += ("--JM", "%dG" % max(1, int(0.8*mem/1024)))
		max_insert = kwargs.get(
							'max_insert',
							diagnostics.lookup_insert_size().max)
		if max_insert:
			self.args += ('--group_pairs_distance', int(max_insert))
		self.add_threading('--CPU')
		self.args += args
		self.run()


class VelvetH (BaseWrapper):
	"""
	velveth component of the Velvet *de novo* assember
	http://www.ebi.ac.uk/~zerbino/velvet/
	"""

	def __init__(self, outdir, kmer, *args, **kwargs):
		self.init('velveth', **kwargs)
		self.version()
		self.max_concurrency = 16
		self.add_openmp()
		self.args += (outdir, kmer) + args
		self.run()


class VelvetG (BaseWrapper):
	"""
	velvetg component of the Velvet *de novo* assember
	http://www.ebi.ac.uk/~zerbino/velvet/
	"""

	def __init__(self, outdir, *args, **kwargs):
		self.init('velvetg', **kwargs)
		self.version()
		self.max_concurrency = 16
		self.add_openmp()
		self.args.append(outdir)
		self.args += args
		self.run()


class VelvetOptimiser (BaseWrapper):
	"""
	Perl script for automatically optimising the three primary parameters
	of the Velvet assembler
	http://bioinformatics.net.au/software.velvetoptimiser.shtml
	"""

	def __init__(self, velveth, *args, **kwargs):
		self.init('VelvetOptimiser.pl', **kwargs)
		self.args.append('-f "%s"' % velveth)
		self.args += args
		self.run()


class Yasra (BaseWrapper):
	"""
	YASRA: comparative assembly of short reads using a reference genome
	http://www.bx.psu.edu/miller_lab
	"""

	def __init__(self, *args, **kwargs):
		self.init('make', **kwargs)
		self.name = 'yasra'
		if not os.path.exists(os.path.join(self.cwd, 'Makefile')):
			utils.die("couldn't find YASRA Makefile in dir '%s'" % self.cwd)
		self.args += args
		self.run()


class YasraMakefile (BaseWrapper):
	"""
	Utility script for generating a Makefile for a YASRA run
	"""

	def __init__(self, reads, template, *args, **kwargs):
		self.init('yasra_makefile', **kwargs)
		self.args += (reads, template)
		self.args += args
		self.run()

# vim: noexpandtab ts=4 sw=4
