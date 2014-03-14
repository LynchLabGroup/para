import parse_gff
import sys, getopt

def usage():
    print sys.argv[1], "-g [--gff] gff file -f [--fasta] fasta file"


def main():

    GFFfile = ''
    fastafile = ''

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hg:f:", ["help", "gff", "fasta"])
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
        
#    print 'GFFfile: ', GFFfile
#    print 'FastaFile: ', fastaFile


    GFFdict = dict() # I'm going to store the parsed GFF lines in a dictionary

    print 'Parsing GFF file ...'
    for record in parse_gff.parseGFF3(GFFfile):
        id = record.attributes["ID"]
        GFFdict[id] = record
        #recordCount += 1
    print 'done.'

    print 'Loading Fasta file in memory ...'
    fastaDict = parse_gff.load_fasta(fastaFile)
    print 'done.'


    #seq = parse_gff.getSequence(GFFdict, fastaDict, "PSEXGNT00001", "CDS")
    #print seq

    prot_seq = parse_gff.get_prot_seq(GFFdict, fastaDict, "PSEXGNT00001")
    print prot_seq



main()

