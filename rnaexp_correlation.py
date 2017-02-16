from genelist_creator import read_file
from scipy.stats.stats import pearsonr
from queue import PriorityQueue
from math import nan

class GeneExpressionSetCorrelation:
    def __init__(self, target_gene, containing_gene, r_row):
        self.target_gene = target_gene
        self.containing_gene = containing_gene
        self.r_row = r_row
    def __lt__(self, other):
        return self.r_row > other.r_row
    def to_string(self):
        return self.target_gene.geneName + " " + self.containing_gene.geneName + " " + str(self.r_row)
        
def main():
    payload = read_file("fullresults.tsv")
    input_gene = payload[0][0]
    find_correlations(input_gene, payload[0])
    
def find_correlations(input_gene, intermediate_gene_list):
    pqueue = PriorityQueue(len(intermediate_gene_list))
    for intermediate_gene in intermediate_gene_list:
        if len(intermediate_gene.rnaSeq) != len(input_gene.rnaSeq):
            continue
        # Determine how to utilize p_val (pearsonr[1])
        r_row = pearsonr(input_gene.rnaSeq, intermediate_gene.rnaSeq)[0]
        # TODO: NAN not replacing
        r_row = r_row if r_row != nan else -1
        pqueue.put((r_row, GeneExpressionSetCorrelation(input_gene, intermediate_gene, r_row)))
    while not pqueue.empty():
        item = pqueue.get()
        # TODO: not printout out in order
        print(item[1].to_string())
        
        
if __name__ == "__main__":
    main()