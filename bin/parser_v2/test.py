import parse_gff_v2
import sys, getopt

def usage():
    print sys.argv[1], "-g [--gff] gff file -f [--fasta] fasta file -l [--list] list of IDs to extract"


def main():

    GFFfile = ''
    fastafile = ''
    listFile = ''

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hg:f:l:", ["help", "gff", "fasta", "list"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-g", "-gff"):
            GFFfile = arg
        elif opt in ("-f", "--fasta"):
            fastaFile = arg
        elif opt in ("-l", "--list"):
            listFile = arg
        
#    print 'GFFfile: ', GFFfile
#    print 'FastaFile: ', fastaFile


#    GFFdict = dict() # I'm going to store the parsed GFF lines in a dictionary

#    a = raw_input()

    print 'Parsing GFF file ...'
    GFFdict = parse_gff_v2.loadGFF(GFFfile, "CDS")
    print 'done.'

#    a = raw_input()

    print 'Loading Fasta file in memory ...'
    fastaDict = parse_gff_v2.loadFasta(fastaFile)
    print 'done.'

#    a = raw_input()

    with open(listFile) as lf:
        for line in lf:
            id = line.rstrip('\n')
            prot_seq = parse_gff_v2.getProteinSeq(GFFdict, fastaDict, id)
            print '>', id
            print prot_seq

#    seq = parse_gff_v2.getSequence(GFFdict, fastaDict, "PSEXGNT00001", "CDS")
#    print seq

#    prot_seq = parse_gff_v2.getProteinSeq(GFFdict, fastaDict, "PSEXGNT00001")
#    print prot_seq



main()

