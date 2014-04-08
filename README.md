# para #
Repository for *Paramecium* Project
This is suppose to be the swissknife for extracting CDSs and finding motifs in them.

## Folders ##

In this repo you can find several repositories, there were created according to insert link:
+ **bin**, standing for **bin**aries, containing all scripts and program,
+ **data**, which is cached in git, it contains the raw data used for the analysis,
+ **doc**, which contains the IPython Notebook that record the project,
+ **results**, it contains all the produced, results produced using data, including figures or analyses.

## bin folder ##

Depending on what you want to do, several programs can be used.

### `bigfoot` folder ###

Contains all the files for parsing bigfoot's outputs

### `gff` folder ###

Contains all the files to extract upstream sequences, CDSs, gene families and deal with GFF files.

### `jobscripts` folder ###

Contains all the jobscripts submitted to Mason. `setup.job` contains the whole process to extract motifs from a given gene family.

### `scripts` folder ###

Folder containing all little scripts to automate certain tasks:

+ **ConvertFastatoPhylip.pl** Perl script to convert Fasta format files to Phylip format -> make sequences usable by PhyML. from [drmuhammadmunir](https://github.com/drmuhammadmunir/perl)
+ **diff.sh** script to compare two (local folders) and print the result in a out file.

   Usage:
   ```shell
   diff.sh dir1 dir2
   ```

+ **editnewick.py** Python script to edit newick trees to fit BigFoot.

   Usage: 
   ```shell
   python editnewick.py inputfile outputfile
   ```

+ **extractmatch.py** Python script to write a fasta file from matching names with another fasta files.

   Usage:
   ```shell
   python extractmatch.py namefile inputfile outputfile
   ```
   
   Where `namefile` is the file from which you want to extract names and `inputfile` the file to match, according to sequence names in `namefile`, extracted sequences from `inputfile` are written in `outputfile`.

+ **fastaheader.py** Python script to edit complex fasta sequences names.
   
   Usage:
   ```shell
   python fastaheader.py fastafile delimiter outputfile
   ```
   
   If your fasta header is:
   ```
   PCAUDG00089|scaffold_0001|151273-155265|-
   ```
   
   And if you want headers with only `PCAUDG00089` you can use the script as follow:
   ```shell
   python fastaheader.py fastafile "|" outputfile
   ```

+ **multialign.sh** creates multialignment file using MUSCLE
+ **mvoverlap.sh** move families described in overlap file
+ **order.sh** recreate sequences file from multialignment to conserve original family order
+ **profiling.py** profile a given program and outputs its profile in a script
+ **tidyup.sh** tidyup directory after use of translatorx
+ **translatorx.sh** use translatorx on a given set of files in a directory
+ **translatorx_vLocal.pl** [TranslatorX](http://translatorx.co.uk/) local version for protein-guided DNA alignment
   
   **NOTICE**: [MUSCLE](http://drive5.com/muscle/) need to be installed (or another multi-alignment program you use) on the machine where you run TranslatorX.

### extract_rand_cds.py ###

This program extract randomly a given number of gene families using a gene family file and retrieves their CDSs, translated or not.
All gene family sequences are put in a single folder, with names according to their family name.

### ntseq.py ###

This program extract all the CDSs sequences of a given gene family size

### stable.py ###

This program correct order in a multialigned file in MUSCLE. Not authored by me see [MUSCLE website](http://drive5.com/muscle/stable.html)

### upstreamlen.py ###

This program load fasta and gff files and write a .csv file of the length of all upstream sequences in all given files.