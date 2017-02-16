from genelist_creator import read_file
from scipy.stats.stats import pearsonr
from queue import PriorityQueue
from math import isnan

class GeneExpressionSetCorrelation:
    def __init__(self, target_gene, containing_gene, r_row, p_val):
        self.target_gene = target_gene
        self.containing_gene = containing_gene
        self.r_row = r_row
        self.p_val = p_val
    def __lt__(self, other):
        return self.r_row < other.r_row
    def to_string(self):
        return self.target_gene.geneName + "," + self.containing_gene.geneName + "," + str(self.r_row)
        
def run_correlation(file, gene):
    payload = read_file(file)
    input_gene = find_gene(gene, payload[0])
    find_correlations(input_gene, payload[0])
    
def find_correlations(input_gene, intermediate_gene_list):
    pqueue = PriorityQueue(len(intermediate_gene_list))
    for intermediate_gene in intermediate_gene_list:
        if len(intermediate_gene.rnaSeq) != len(input_gene.rnaSeq):
            continue
        # Determine how to utilize p_val (pearsonr[1])
        r_row, p_val = pearsonr(input_gene.rnaSeq, intermediate_gene.rnaSeq)
        r_row = r_row if not isnan(r_row) else -1
        pqueue.put((-r_row, GeneExpressionSetCorrelation(input_gene, intermediate_gene, r_row, p_val)))
    count = 0
    while not pqueue.empty():
        item = pqueue.get()
        if item[1].r_row < 0.8:
            break
        print(item[1].to_string())
        count += 1
def find_gene(gene_name, gene_list):
    for gene in gene_list:
        if gene.geneName == gene_name:
            return gene
        