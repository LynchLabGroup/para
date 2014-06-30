# para #
Repository for *Paramecium* Project
This is suppose to be the swissknife for extracting CDSs and finding motifs in them.
To understand how to run the pipeline, please read the `HOWTO.md` file, the `README.md` file describes the architecture and the structure of the project.

## Folders ##

In this repo you can find several repositories, there were created according to insert link:
+ **bin**, standing for **bin**aries, containing all scripts and program,
+ **data**, which is cached in git, it contains the raw data used for the analysis,
+ **doc**, which contains the IPython Notebook that record the project,
+ **results**, it contains all the produced, results produced using data, including figures or analyses.

## bin folder ##

Depending on what you want to do, several programs can be used.

### `bigfoot` folder ###

Contains all the files for parsing bigfoot's outputs:

+ **classify.py** command-line interface script to produce various indexes on motifs
   Usage:
   ```shell
   python classify.py seqtabfile outputfile -f seqfieldnumber [--header]
   ```
   see `python classify.py -h` for more detailed help

+ **memecomp.py** use Biopython motifs subpackage to compute the differences between MEME and BigFoot outputs

   Usage:
   ```shell
   python memecomp.py bigfoot_output meme_output
   ```
   See `python memecomp.py -h` for detailed help

+ **memeextract.py** extract motifs from MEME XML files.
   Usage:
   ```shell
   python memeextract.py evalue outputfile memexmlfile1 [memexmlfile2 ...]
   ```
   See `python memeextract.py -h` for detailed help

+ **parser.py** contains the class SeqParser to parse BigFoot's outputs file `.mpd` and `.pred`

+ **setup.py** uses `weight.py` and `parser.py` to parse BigFoot's outputs file. Command-line interface.
   Usage:
   ```shell
   python setup.py bigfoot_file.mpd bigfoot_file.pred
   ```
   See `python setup.py -h` for detailed help
+ **stamp.py** generate motif alignment from `memecomp.sh` output
   Usage:
   ```shell
   python stamp.py memecomp_output output_file
   ```
+ **weight.py** contains WeightedSeq class, to parse the various sequences in `.pred` files

### `gff` folder ###

Contains all the files to extract upstream sequences, CDSs, gene families and deal with GFF files.

### `jobscripts` folder ###

Contains all the jobscripts submitted to Mason. `setup.job` contains the whole process to extract motifs from a given gene family.

### `scripts` folder ###

Folder containing all little scripts to automate certain tasks:

+ **allcomp.sh** contains command to create a summary of found motifs in .comp files

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

+ **erasediscard.sh** erase families in discarded files, *USELESS* now that discarded families are not retrieved

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
   ``` motifs
   
   And if you want headers with only `PCAUDG00089` you can use the script as follow:
   ```shell
   python fastaheader.py fastafile "|" outputfile
   ```

+ **genomeextract.sh** use custom script to extract motifs in all species, work on BigFoot's and MEME's

+ **highnonmatch.R**, R script to extract non ribosomal highly expressed genes 

+ **intergenic.sh** script to compute intergenic distances in specific species
   **/!\\** Modify to change path to files.

+ **intmotifs.sh** script to extract motifs sizes and distances from results.
+ **memecomp.sh** script to use MEME and `memecomp.py` on all available families
+ **memelen.sh** show the length of all motifs in `.meme.motifs` files
+ **motifslengths.sh** command example to compute the length of all BigFoot's motifs
+ **multialign.sh** creates multialignment file using MUSCLE
+ **mvoverlap.sh** move families described in overlap file
+ **order.sh** recreate sequences file from multialignment to conserve original family order
+ **profiling.py** profile a given program and outputs its profile in a script
+ **swaporder.py** swap sequence order in fasta file

   Usage:
   ```shell
   python swaporder.py input_fasta output_fasta
   ```
   
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

## results folder ##

Contains data processed from raw data. `analysis` contains scripts and output files of analysis.

### analysis folder ###

Contains several scripts to study various subjects:

+ **familysize** contains scripts to study family sizes and according results
+ **intergene** contains intergenic distance computation
+ **motifs sizes** contains everything to compute sizes of motifs and family

### TO DO ###

- [ ] Implement testing
- [ ] Add example codes
- [ ] Tidy Up repo (delete tests and unused bash script)
- [X] Implement R script CLI
- [X] Implement bash script CLI 