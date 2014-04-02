# para #
Repository for *Paramecium* Project
This is suppose to be the swissknife for extracting CDSs and finding motifs in them.

## Folders ##

In this repo you can find several repositories, there were created according to insert link:
+ **bin**, standing for **bin**aries, containing all scripts and program,
+ **data**, which is cached in git, it contains the raw data used for the analysis,
+ **doc**, which contains the IPython Notebook that record the project,
+ **results**, it contains all the produced, results produced using data, including figures or analyses.

## Programs ##

Depending on what you want to do, several programs can be used.

### extract_rand_cds ###

This program extract randomly a given number of gene families using a gene family file and retrieves their CDSs, translated or not.
All gene family sequences are put in a single folder, with names according to their family name.

### ntseq ###

This program extract all the CDSs sequences of a given gene family size

### 