# Program to output a file
import parse_gff_v2
import sys, operator

def main():

    gffFile = sys.argv[1]
    fastaFile = sys.argv[2]
    spec = sys.argv[3]


    sys.stderr.write("Load fasta file in memory ...")
    fastaDict = parse_gff_v2.load_fasta(fastaFile)
    sys.stderr.write("done.\n")

    sys.stderr.write("Parsing GFF file ...")
    GFFdict = parse_gff_v2.loadGFFbySeqID(gffFile, "mRNA")
    sys.stderr.write("done.\n")


    for scafID in GFFdict.keys():
        GFFdict[scafID].sort(key=operator.attrgetter("start"))
        nbEntries = len(GFFdict[scafID])

        for i in range(0, nbEntries):
            p = True
            intergenicSize = -1
            currentGene = GFFdict[scafID][i]
            if i > 0 and i < (nbEntries-1):
                previousGene = GFFdict[scafID][i-1]
                intergenicSize = (currentGene.start - previousGene.end) + 1
            elif i == 0:
                intergenicSize = currentGene.start
            elif i == (nbEntries-1):
                scafSize = len(fastaDict[scafID])
                intergenicSize = (scafSize - currentGene.end)
                if intergenicSize == 0:
                    p = False
            if p:
                print currentGene.attributes["Name"], intergenicSize, spec


if __name__ == "__main__":
    main()
