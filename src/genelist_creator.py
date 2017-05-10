import csv
from pprint import pprint

class GeneRnaDataSet:
    def __init__(self, dbIdentifier='', rnaSeq = []):
        self.dbIdentifier = dbIdentifier
        self.rnaSeq = rnaSeq
    def __str__(self):
        return self.dbIdentifier + " " +  str(self.rnaSeq)

class GeneIdentifierSet:
    def __init__(self, dbIdentifier='', secondaryIdentifier='', geneName='', synonyms = []):
        self.dbIdentifier = dbIdentifier
        self.secondaryIdentifier = secondaryIdentifier
        self.geneName = geneName
        self.synonyms = synonyms
        
    
# Works    
# Unmarshalls a TSV file containing gene's expression stage values
#
# Takes in a TSV file with 5 columns
# COL 0 : Gene DB Identifier
# COL 1: Gene Secondar Identifier
# COL 2: Gene Name
# COL 3: RNA Expression Stage
# COL 4: RNA Expression Score
#
# Takes in a TSV file with 3 columns
# COL 0 : Gene DB Identifier
# COL 1: RNA Expression Stage
# COL 2: RNA Expression Score
#
# Returns a List of GeneExpressionSet objects 
#
# This method usually acts as a subroutine for another method, not a standolone function
def read_file_rnadata(filename):
    with open(filename,'rt') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        gene = GeneRnaDataSet()
        geneList = []
        rnaSeq = []
        count = 0
        geneString = ""
        for row in tsvin:
            if geneString != row[0]:
                gene.rnaSeq = rnaSeq
                geneList.append(gene)
                gene = GeneRnaDataSet()
                rnaSeq = []
                gene.dbIdentifier = geneString = row[0]
                
            rnaSeq.append(int(row[2]))
            count += 1
        geneList = geneList[1:]
        return [geneList, 0, max]


# Unmarshalls a TSV file containing gene's expression stage values
#
# Takes in a TSV file with 4 columns
# COL 0 : Gene DB Identifier
# COL 1: Gene Secondar Identifier
# COL 2: Gene Name
# COL 3: Gene Synonym
#
# Returns a List of GeneExpressionSet objects 
#
# This method usually acts as a subroutine for another method, not a standolone function
def read_file_geneid(filename):
    with open(filename,'rt') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        gene = GeneIdentifierSet()
        geneList = []
        synonyms = []
        count = 0
        geneString = ""
        for row in tsvin:
            if geneString != row[0]:
                gene.synonyms = synonyms
                geneList.append(gene)
                gene = GeneIdentifierSet()
                synonyms = []
                gene.dbIdentifier = geneString = row[0]
                gene.secondaryIdentifier = row[1]
                gene.geneName = row[2]
                
            synonyms.append(int(row[4]))
            count += 1
        geneList = geneList[1:]
        return [geneList, 0, max]

if __name__ == '__main__':
    # This file contains all RNA expression data for all genes available in D. m
    #read_file("../res/results_w_more_ids.tsv") To access deprecated read_file method, go to commit #db4a6ee
    
    for geneRnaDataSet in read_file_rnadata("../res/tsv_files/rna_data_limited_ids.tsv"):
        print(geneRnaDataSet.dbIdentifier)
    #pprint(read_file_geneid("../res/tsv_files/gene_identifiers.tsv"))          
    
    # This will not work on new code since it lacks all required columns. Refer to old version of code
    # read_file("../res/results_w_more_ids.tsv")
    
    
    
    
    
    
    
    
    
    
    