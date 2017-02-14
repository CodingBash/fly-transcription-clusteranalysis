import csv

import math
import random

# Only needed if you want to display your plots inline if using Notebook
# change inline to auto if you have Spyder installed

plotly = False
try:
    import plotly
    from plotly.graph_objs import Scatter, Scatter3d, Layout
except ImportError:
    print("INFO: Plotly is not installed, plots will not be generated.")


def main():
    fileResults = readFile();
    lower = fileResults[1]
    upper = fileResults[2]
    

class GeneExpressionSet:
    def __init__(self):
        self.geneName = ''
        self.rnaSeq = []

def readFile():
    with open('fullresults.tsv','rt') as tsvin:
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
        max = 0
        geneMax = GeneExpressionSet()
        for gene in geneList:
            if len(gene.rnaSeq) != 104:
                print(gene.geneName)
                print(len(gene.rnaSeq))
            for rnaSeq in gene.rnaSeq:
               if rnaSeq > max:
                   max = rnaSeq
                   geneMax = gene
        print("MAX " + str(max))
        print(geneMax.geneName)
        print(geneMax.rnaSeq)
        print(len(geneList))
        print(count)
        # Interesting... dimensions either 74 0r 104
        return [geneList, 0, max]


def cluster(geneList):
    print("Hello World")

if __name__ == "__main__":
    main()