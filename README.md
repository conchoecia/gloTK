# gloTK

GloTK is the Genomes of Luminous Organisms Toolkit. It contains
scripts to facilitate genome assembly optomization and analysis using
common genome assemblers.

## Supported Assemblers

Currently there are only scripts for one assembler:

- Meraculous *2.2.2.4+*

## Dependencies

- The `mer_reporter.py` script currently requires the following to be installed:
  - Python 3.5
    - (If you haven't installed yet, I recommend the [Anaconda](https://www.continuum.io/downloads) distribution)
  - [py-gfm](https://py-gfm.readthedocs.io/en/latest/)
    - _Python Github-flavored Markdown_

## Install

You can install the `gloTK` package on your linux machine by typing
`pip install gloTK` in your terminal.

Major releases of `gloTK` are uploaded to the [Python Package Index](https://pypi.python.org/pypi/gloTK).

## Scripts

To use these programs, install `gloTK` as described above and type `glotk` into
your terminal. If you hit the `tab` key a few times your computer should try to
autocomplete and you will see a few options:

- **glotk-mer-reporter.py**
  - Generates HTML reports of all Meraculous assemblies found within
    the current working directory _(CWD)_. Collates all relevant files and images
    into a `meraculous_reports` directory in the _CWD_.
- **glotk-sweep**
  - This allows the user to input a single meraculous config file, and to define 
    parameters to sweep on and optomize. Currently this program supports sweeping 
    on `mer_size` and `bubble_depth_threshold`. Typing `glotk-sweep` into your
    terminal after installing will give you more details.

## News


### _*v0.1.15*_ - 20161116


* This version incorporates glotk-project. This script reads in a config file
  and generates a 'glotk_project' directory in which information about the
  libraries, assemblies, et cetera for this Meraculous genome project is stored. 
  Assemblies, reports, config files, et cetera  will take place in
  adjacent directories
* gloTK project files will now take on the following file structure:

```
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
|--glotk_assemblies
|  |--as000_<assembly_name>/
|  |  |--assembler_output/
|  |
|  |--as001_20161102_ME_pleu_bach_k63_d1/
|     |--meraculous_etcetera
|
|--glotk_configs
|     as000_<assembly_name>.config
|     as001_20161102_ME_pleu_bach_k63_d1.config
|
|--glotk_reads
|  |  reads.log
|  |
|  |--reads0/
|  |     <forward_reads_symlink>.fq.gz
|  |     <reverse_reads_symlink>.fq.gz
|  |
|  |--reads1/
|        <processed_forward>.fq.gz
|        <processed_reverse>.fq.gz
|
|--glotk_fastqc
|--glotk_kmers
|--glotk_reports
```


### _*v0.1.14*_ - 20161102


* Gave a --triplet flag to glotk-sweep to cause the assembler to run diploid modes
  0, 1, and 2 on the same kmer size.


### _*v0.1.13*_ - 20161102

* Made glotk-sweep function by inputting each value of k to assemble for, rather
  than using 'sweep start', 'sweep stop', and 'sweep interval'.
  Usage is: '--slist 23 27 57' to perform assemblies for those three values
