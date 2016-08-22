#!/usr/bin/env python3

# Tutorials on how to setup python package here:
#   - http://python-packaging.readthedocs.io/en/latest/testing.html
#   - https://jeffknupp.com/blog/2013/08/16/open-sourcing-a-python-project-the-right-way/

from setuptools import setup

setup(name='gloTK',
      version='0.1',
      description='Genomes of Luminous Organisms Toolkit',
      long_description="""
          GloTK is a toolkit for enhancing productivity in the process
          of whole genome de novo assembly. Currently the package
          supports automated assembly HTML report generation and automated
          assemly parameter optomization for the Meraculous assembler.""",
      url='https://github.com/cypridina/gloTK',
      author='Darrin Schultz',
      author_email='dschultz@mbari.org',
      classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Operating System :: POSIX :: Linux',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Intended Audience :: Science/Research'
          ],
      license='GPLv3',
      provides=['gloTK'],
      packages=['gloTK', 'gloTK.scripts'],
      install_requires=[
          #MerParse - None
          #tests - nose (see below)
          #MerRunAnalyzer
          "py_gfm",
          "pymdown-extensions",
          "markdown"
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      entry_points = {
            'console_scripts': ['glotk-sweep=gloTK.scripts.glotk_sweep:main'],
        },
      zip_safe=False,
      include_package_data=True)
