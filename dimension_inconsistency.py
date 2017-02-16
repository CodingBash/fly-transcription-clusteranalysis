from genelist_creator import read_file, GeneExpressionSet



def main():
    payload = read_file('fullresults.tsv')
    geneList = payload[0]
    max = 0
    geneMax = GeneExpressionSet()
    lowDimensionalGeneList = []
    for gene in geneList:
        if len(gene.rnaSeq) != 104:
            lowDimensionalGeneList.append(gene)
        for rnaSeq in gene.rnaSeq:
           if rnaSeq > max:
               max = rnaSeq
               geneMax = gene

    for gene in lowDimensionalGeneList:
        print(gene.geneName)
    print("COUNT: " + str(len(lowDimensionalGeneList)))
    print("MAX " + str(max))
    print(geneMax.rnaSeq)
    print(geneMax.geneName)

if __name__ == "__main__":
    main()