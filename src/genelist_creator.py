import csv

class GeneExpressionSet:
    def __init__(self, dbIdentifier='', secondaryIdentifier='', geneName='', rnaSeq = []):
        self.dbIdentifier = dbIdentifier
        self.secondaryIdentifier = secondaryIdentifier
        self.geneName = geneName
        self.rnaSeq = rnaSeq

# Unmarshalls a TSV file containing gene's expression stage values
#
# Takes in a TSV file with 5 columns
# COL 0 : Gene DB Identifier
# COL 1: Gene Secondar Identifier
# COL 2: Gene Name
# COL 3: RNA Expression Stage
# COL 4: RNA Expression Score
#
# Returns a List of GeneExpressionSet objects 
#
# This method usually acts as a subroutine for another method, not a standolone function
def read_file(filename):
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
                gene.dbIdentifier = geneString = row[0]
                gene.secondaryIdentifier = row[1]
                gene.geneName = row[2]
                
            rnaSeq.append(int(row[4]))
            count += 1
        geneList = geneList[1:]
        return [geneList, 0, max]

if __name__ == '__main__':
    # This file contains all RNA expression data for all genes available in D. m
    read_file("../res/results_w_more_ids.tsv")
    
    # This will not work on new code since it lacks all required columns. Refer to old version of code
    # read_file("../res/results_w_more_ids.tsv")