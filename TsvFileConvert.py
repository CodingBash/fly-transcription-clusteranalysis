import csv

class GeneExpressionSet:
    def __init__(self):
        self.geneName = ''
        self.rnaSeq = []

with open('fullresults.tsv','rt') as tsvin:
    # Read TSV file
    tsvin = csv.reader(tsvin, delimiter='\t')

    # Setup initial loop variables
    gene = GeneExpressionSet()
    geneList = []
    rnaSeq = []
    count = 0

    # Incrememt thru each row in TSV
    for row in tsvin:
        # If a new gene is encountered in row
        if count % 104 == 0:
           # Add gene to list
           gene.rnaSeq = rnaSeq
           geneList.append(gene)
           # Initialize a new gene
           gene = GeneExpressionSet()
           rnaSeq = []
           gene.geneName = row[0]
        rnaSeq.append(row[2])
        count += 1
    # Trim first entry of geneList (since empty)
    geneList = geneList[1:]

    
    for gene in geneList:
        print(gene.geneName)
        print(gene.rnaSeq)

    # Print data
    print(len(geneList))