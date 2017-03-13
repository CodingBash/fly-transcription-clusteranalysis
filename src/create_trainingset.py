import json
import csv
from pprint import pprint

tf_dict = dict()
with open("../res/training_set_data/tf_to_genes_results.tsv", 'rt') as tf_gene_tsvin:
    tf_gene_tsvin = csv.reader(tf_gene_tsvin, delimiter='\t')
    for row in tf_gene_tsvin:
        tf = row[1]
        gene = row[3]
        if not tf_dict.get(tf):
            tf_dict.update({tf: list()})
        tf_list = tf_dict.get(tf)
        tf_list.append(gene)
        tf_dict.update({tf: tf_list})

gene_dict = dict()
with open("../res/training_set_data/gene_to_tf_results.tsv", 'rt') as gene_tf_tsvin:
    gene_tf_tsvin = csv.reader(gene_tf_tsvin, delimiter='\t')
    for row in gene_tf_tsvin:
        gene = row[1]
        tf = row[3]
        if not gene_dict.get(gene):
            gene_dict.update({gene: dict()})
        gene_tf_dict = gene_dict.get(gene)
        gene_tf_dict.update({tf: tf_dict.get(tf)})
        gene_dict.update({gene: gene_tf_dict})

with open("../res/training_set_data/training_set_results.json", 'w') as jsonfile:
    json.dump(gene_dict, jsonfile, indent=5)
    
pprint("JSON contains " + str(len(gene_dict)) + " genes.")