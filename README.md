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
