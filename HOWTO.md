# How to run the pipeline? #

## Pipeline Requirements ##

Be sure everything is installed and accessible

- Python 2.7
- [Biopython (>=1.63)](http://biopython.org/)
- Java (compatible with BigFoot use)
- [BigFoot](http://sourceforge.net/projects/bigfoot/) (in the project's parent folder)
- Perl
- [TranslatorX](http://translatorx.co.uk/) local version (in the `scripts`folder)
- [Muscle](http://www.drive5.com/muscle/) (to use within TranslatorX, change TranslatorX options to use another alignment program)
- [PhyML](http://www.atgc-montpellier.fr/phyml/binaries.php)

## What are the files to run the pipeline? ##

The pipeline consists of a `Makefile` and a configuration file `config.mak`.
To have more informations on Make, see the [official documentation](https://www.gnu.org/software/make/manual/).

### `Makefile` ###
This file contains the proper pipeline with rules to run it.
**DO NOT MODIFY** it, unless you know what you are doing.

### `config.mak` ###
This file contains all the parameters, **modify `config.mak`** for particular uses of the pipeline.

## What commands to use? ##

**/!\\** All the genomes and GFF annotation files are hard-coded in Python files, if you want to use other genomes/annotations, modify `bin/ntseq.py` and `bin/gff/main.py` **/!\\**

**First**, we want to retrieve all the upstream sequences and their according CDSs.
Modify the `config.mak` to adjust parameters, upstream sequences will be put in the `UP` folder specified in `config.mak` and CDSs in `CDS` folder of `config.mak`. Unsure these folders are created before using the pipeline.

Run while in the "para/" directory
```shell
make retrieve_upstream
make retrieve_CDS
```
These commands won't run on parallel to only load genomes sequences and annotations once.

**Then**, to launch the pipeline run:
```shell
make all -r
```
the `-r` option tells make not to use default rules, as this Makefile is twisted from the original uses of `make`.

When submitting a job using TORQUE, you can parallelize the pipeline by using the `-j` flag.
Example of Jobscript:
```shell
#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=4,vmem=18gb,walltime=02:00:00:00
#PBS -M john.doe@example.com
#PBS -m be
#PBS -N JobName
#PBS -j oe
make all -r -j 4
```

### How do I run a single step? ###

If you want to generate specific files, you can call them directly:
```shell
make results/WGD2ANC00002/WGD2ANC00002.comp -r
```

This will run the approriate rules to make `WGD2ANC00002.comp`

For the moment, it is not possible to rerun simply a given step in the pipeline without specifying all the file names.
