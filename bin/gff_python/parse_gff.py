#!/usr/bin/env python
"""
A simple parser for the GFF3 format.

Test with transcripts.gff3 from
http://www.broadinstitute.org/annotation/gebo/help/gff3.html.

Format specification source:
http://www.sequenceontology.org/gff3.shtml

Version 1.0
"""
from __future__ import with_statement
from collections import namedtuple
import gzip
import urllib

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC



__author__  = "Uli Koehler"
__license__ = "Apache License v2.0"

#Initialized GeneInfo named tuple. Note: namedtuple is immutable
gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)

def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dict"""#
    if attributeString == ".": return {}
    ret = {}
    for attribute in attributeString.split(";"):
        key, value = attribute.split("=")
        ret[urllib.unquote(key)] = urllib.unquote(value)
    return ret

def parseGFF3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.
    
    Supports transparent gzip decompression.
    """

    #Parse with transparent decompression
    openFunc = gzip.open if filename.endswith(".gz") else open
    with openFunc(filename) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalizedInfo = {
                "seqid": None if parts[0] == "." else urllib.unquote(parts[0]),
                "source": None if parts[1] == "." else urllib.unquote(parts[1]),
                "type": None if parts[2] == "." else urllib.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urllib.unquote(parts[6]),
                "phase": None if parts[7] == "." else urllib.unquote(parts[7]),
                "attributes": parseGFFAttributes(parts[8])
            }
            #Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            yield GFFRecord(**normalizedInfo)


#######################################################################################################################################################
#
#  JF's part starts here

            

# getSequence : return the DNA sequence of a gene/transcript/mRNA based on its parent ID
# Arguments:
#  - GFFdict: a dictionary containing all the GFF lines parsed into objects (key=seqid / value=object returned by the GFF parser for one line of GFF)
#  - fastaDict: a simple dictionary of the fasta file (key=sequence ID / value=sequence)
#  - id: the parent ID (typically: PSEXGNT00001)
#  - type: the sequence type from GFF (typically one of 'CDS', 'exon', ...
#
# Returns the DNA sequence
#
# WARNING: does not do any asserts. Passing an id that does not have any parent will crash the programm...
def get_seq(GFFdict, fastaDict, seq_id, seq_type):
    """return the DNA sequence of a gene/transcript/mRNA based on its parent ID."""
    ll = list()

    for key, gffEntry in GFFdict.items():
        if gffEntry.type == seq_type:
            if gffEntry.attributes["Parent"] == seq_id:
                ll.append(gffEntry)

    ll.sort(key=lambda x: x.start, reverse=False)

    seq = ""
    for entry in ll:        
        exon_seq = fastaDict[entry.seqid][(entry.start-1):(entry.end)]
        seq += exon_seq

    bseq = Seq(str(seq), Bio.Alphabet.generic_dna)
    if entry.strand == "-":
        bseq = bseq.reverse_complement()

    seq = str(bseq)
    return seq


#####################################################################################################################################################
# getProteinSeq : return the translated sequence of a gene/transcript/mRNA based on its parent ID
# Arguments:
#  - GFFdict: a dictionary containing all the GFF lines parsed into objects (key=seqid / value=object returned by the GFF parser for one line of GFF)
#  - fastaDict: a simple dictionary of the fasta file (key=sequence ID / value=sequence)
#  - id: the parent ID (typically: PSEXGNT00001)
#
# Returns the translated (protein) sequence
#
# WARNING: does not do any asserts. Passing an id that does not have any parent will crash the programm...
def get_prot_seq(GFFdict, fastaDict, id, table=6):
    """Return the translated sequence of a gene/transcript/mRNA based on its parent ID"""
    dna_seq = Seq(get_seq(GFFdict, fastaDict, id, "CDS"), Bio.Alphabet.generic_dna)
    prot_seq = dna_seq.translate(table=table)
    return str(prot_seq)
    

######################################################################################
# loadFasta: load a fasta file into a dictionary (key=sequence ID / value=sequence)
# Arguments:
# - in_file: path to the fasta file
#
# WARNING: does not do any check on the fasta file (providing a non-fasta file may crash the programm... or even worse, let it go without warning and lead to troubles downstream)
def load_fasta(in_file):
    """Return a dictionnary with sequences indexed by sequence id."""
    fastaDict = dict()
    with open(in_file, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta") :
            fastaDict[record.id] = record.seq
    return fastaDict
        


