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
    def set_rank(self, rank):
        self.rank = rank
    def to_string(self):
        return self.target_gene.geneName + "," + self.containing_gene.geneName + "," + str(self.r_row)
        
def run_correlation(file, gene, corr_limit=0.3, output_file=""):
    payload = read_file(file)
    input_gene = find_gene(gene, payload[0])
    find_correlations(input_gene, payload[0], corr_limit, output_file)
    
            
def run_sanitized_correlation(file, gene, corr_limit=0.3, output_file=""):
    payload = read_file(file)
    input_gene = find_gene(gene, payload[0])
    find_sanitized_correlations(input_gene, payload[0], corr_limit, output_file)
   
# TODO: Fix the risk of NoneType
# TODO: Add rank    
def find_correlations(input_gene, intermediate_gene_list, corr_limit, output_file=""):
    output_to_file = len(output_file) > 0
    pqueue = PriorityQueue(len(intermediate_gene_list))
    for intermediate_gene in intermediate_gene_list:
        if len(intermediate_gene.rnaSeq) != len(input_gene.rnaSeq):
            continue
        r_row, p_val = pearsonr(input_gene.rnaSeq, intermediate_gene.rnaSeq)
        r_row = r_row if not isnan(r_row) else -1
        pqueue.put((-r_row, GeneExpressionSetCorrelation(input_gene, intermediate_gene, r_row, p_val)))
    count = 0
    if output_to_file:
        n_file = open(output_file, 'w')
    while not pqueue.empty():
        item = pqueue.get()
        if item[1].r_row < corr_limit:
            break
        if output_to_file:
            n_file.write(item[1].to_string() + "\n")
        else:
            print(item[1].to_string())
        count += 1
    if output_to_file:
        n_file.close();
                    
# TODO: Fix the risk of NoneType
# TODO: Add rank    
def find_sanitized_correlations(input_gene, intermediate_gene_list, corr_limit, output_file=""):
    output_to_file = len(output_file) > 0
    pqueue = PriorityQueue(len(intermediate_gene_list))
    
    #Sanitize input_gene rnaSeq data
    input_gene_rnaSeq = []
    for index in range(0, len(input_gene.rnaSeq)):
        if input_gene.rnaSeq[index] != 0:
            input_gene_rnaSeq.append([index, input_gene.rnaSeq[index]])
                
    for intermediate_gene in intermediate_gene_list:
        if len(intermediate_gene.rnaSeq) != len(input_gene.rnaSeq):
            continue
        
        # From input_gene sanitized data, sanitize the intermediate gene
        input_gene_rnaSeq_index = 0
        intermediate_gene_rnaSeq = []
        for index in range(0, len(intermediate_gene.rnaSeq)):
            if index == input_gene_rnaSeq[input_gene_rnaSeq_index][0]:
                intermediate_gene_rnaSeq.append([index, intermediate_gene.rnaSeq[index]])
                input_gene_rnaSeq_index += 1
        
        # Take pearsonr of tuples' second column (1st when 0-based)
        r_row, p_val = pearsonr([item[1] for item in input_gene_rnaSeq], [item[1] for item in intermediate_gene_rnaSeq])
        r_row = r_row if not isnan(r_row) else -1
        pqueue.put((-r_row, GeneExpressionSetCorrelation(input_gene, intermediate_gene, r_row, p_val)))
    count = 0
    if output_to_file:
        n_file = open(output_file, 'w')
    while not pqueue.empty():
        item = pqueue.get()
        if item[1].r_row < corr_limit:
            break
        if output_to_file:
            n_file.write(item[1].to_string() + "\n")
        else:
            print(item[1].to_string())
        count += 1
    if output_to_file:
        n_file.close();

# TODO: Fix the risk of NoneType
def find_correlations_simple(input_gene, intermediate_gene_list, corr_limit):
    pqueue = PriorityQueue(len(intermediate_gene_list))
    for intermediate_gene in intermediate_gene_list:
        if len(intermediate_gene.rnaSeq) \
              != len(input_gene.rnaSeq):
            continue
        r_row, p_val = pearsonr(input_gene.rnaSeq, intermediate_gene.rnaSeq)
        r_row = r_row if not isnan(r_row) else -1
        pqueue.put((-r_row, GeneExpressionSetCorrelation(input_gene, intermediate_gene, r_row, p_val)))
    
    correlation_list = []
    count = 0
    while not pqueue.empty():
        item = pqueue.get()
        if item[1].r_row < corr_limit:
            break
        correlation_list.append([item[1], count])
        count += 1
    return correlation_list

# TODO: Fix the risk of NoneType
def find_sanitized_correlations_simple(input_gene, intermediate_gene_list, corr_limit):
    pqueue = PriorityQueue(len(intermediate_gene_list))
    
    #Sanitize input_gene rnaSeq data
    input_gene_rnaSeq = []
    for index in range(0, len(input_gene.rnaSeq)):
        if input_gene.rnaSeq[index] != 0:
            input_gene_rnaSeq.append([index, input_gene.rnaSeq[index]])
                
            
            
    for intermediate_gene in intermediate_gene_list:
        if len(intermediate_gene.rnaSeq) \
              != len(input_gene.rnaSeq):
            continue
        
        # From input_gene sanitized data, sanitize the intermediate gene
        input_gene_rnaSeq_index = 0
        intermediate_gene_rnaSeq = []
        for index in range(0, len(intermediate_gene.rnaSeq)):
            if index == input_gene_rnaSeq[input_gene_rnaSeq_index][0]:
                intermediate_gene_rnaSeq.append([index, intermediate_gene.rnaSeq[index]])
                if input_gene_rnaSeq[input_gene_rnaSeq_index][0] != input_gene_rnaSeq[len(input_gene_rnaSeq)-1][0]:
                    input_gene_rnaSeq_index += 1
                
                
        # Take pearsonr of tuples' second column (1st when 0-based)
        r_row, p_val = pearsonr([item[1] for item in input_gene_rnaSeq], [item[1] for item in intermediate_gene_rnaSeq])
        r_row = r_row if not isnan(r_row) else -1
        pqueue.put((-r_row, GeneExpressionSetCorrelation(input_gene, intermediate_gene, r_row, p_val)))
    
    correlation_list = []
    count = 0
    while not pqueue.empty():
        item = pqueue.get()
        if item[1].r_row < corr_limit:
            break
        correlation_list.append([item[1], count])
        count += 1
    return correlation_list
        
def find_gene(gene_name, gene_list):
    for gene in gene_list:
        if gene.geneName == gene_name:
            return gene
        