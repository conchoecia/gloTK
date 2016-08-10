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

You can execute the scripts by their names if you add the install directory to your `$PATH` variable. For example, in my `~/.bash_profile` file, I added this line to the end:

```
PATH=$PATH:/home/<my_username>/git/gloTK/scripts
```

## Scripts

- *mer_reporter.py*
  - Generates HTML reports of all Meraculous assemblies found within
    the current working directory _(CWD)_. Collates all relevant files and images
    into a `meraculous_reports` directory in the _CWD_.
