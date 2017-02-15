import csv

class GeneExpressionSet:
    def __init__(self):
        self.geneName = ''
        self.rnaSeq = []

def readFile(filename):
    with open(filename,'rt') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        gene = GeneExpressionSet()
        geneList = []
        rnaSeq = []
        count = 0
        geneString = ""
        for row in tsvin:
            if geneString != row[0]:
                gene.rnaSeq = rnaSeq
                geneList.append(gene)
                gene = GeneExpressionSet()
                rnaSeq = []
                gene.geneName = geneString = row[0]
                
            rnaSeq.append(int(row[2]))
            count += 1

        geneList = geneList[1:]
        print(len(geneList))
        print(count)
        # Interesting... dimensions either 74 0r 104
        return [geneList, 0, max]
