from rnaexp_correlation import find_correlations_simple, find_gene, GeneExpressionSetCorrelation
from genelist_creator import read_file, GeneExpressionSet


def run_group_correlation(file, gene_list, corr_limit=0.3, output_file=""):
    output_to_file = len(output_file) > 0
    payload = read_file(file)
    gene_correlation_group = []
    for i in range(len(gene_list)):
        print(gene_list[i])
        input_gene = find_gene(gene_list[i], payload[0])
        if input_gene != None:
            correlation_list = find_correlations_simple(input_gene, payload[0], corr_limit)
            for k in range(len(gene_list)):
                if i == k:
                    continue
                target_gene = find_gene(gene_list[k], payload[0])
                if target_gene != None:
                    found = False
                    for correlated_entry in correlation_list:
                        if correlated_entry[0].containing_gene == target_gene:
                            gene_correlation_group.append([GeneExpressionSetCorrelation(input_gene, target_gene, correlated_entry[0].r_row, correlated_entry[0].p_val), correlated_entry[1]])
                            found = True
                    if found == False:
                        gene_correlation_group.append([GeneExpressionSetCorrelation(input_gene, target_gene, -1, 0), -1])
                else:
                    gene_correlation_group.append([GeneExpressionSetCorrelation(input_gene, GeneExpressionSet(gene_list[k]), -1, 0), -1])
        else:
            gene_correlation_group.append([GeneExpressionSetCorrelation(GeneExpressionSet(gene_list[i]), GeneExpressionSet(gene_list[i]), -1, 0), -1])

    if output_to_file:
        n_file = open(output_file, 'w')
        
    for gene_correlation in gene_correlation_group:
        if output_to_file:
            n_file.write(gene_correlation[0].to_string() + "," + str(gene_correlation[1]) + '\n')
        else:
            print(gene_correlation.to_string())