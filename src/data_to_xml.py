from genelist_creator import read_file
from xml.etree.ElementTree import Element, SubElement, Comment, tostring

# Creates XML from data where each rna stage is a new child element
def create_xml(input_filename, output_filename, is_prettify=True, gene_count=-1 ):
    payload = read_file(input_filename)
    data = payload[0]
    root = Element("genome")
    root.set('version', '1.0.0')
    root.append(Comment("Generated XML"))
    outer_count = 0;
    for entry in data:
        if gene_count != -1 and outer_count > gene_count:
            break
        if len(entry.rnaSeq) == 104:
            gene_child = SubElement(root, "gene")
            gene_name_child = SubElement(gene_child, "gene-name")
            gene_name_child.text = entry.geneName
            rna_exp_child = SubElement(gene_child, "rna-exp")
            count = 0;
            for rna_exp_stage in entry.rnaSeq:
                rna_exp_stage_child = SubElement(rna_exp_child, "rna-exp-stage", {"order": str(count) })
                rna_exp_stage_child.text = str(rna_exp_stage)
                count+=1
            outer_count += 1
    file = open(output_filename, 'w')
    if is_prettify:
        file.write(prettify(root))
    else:
        file.write(str(tostring(root)))
    file.close()

# Creates XML from data where each rna stage is in the same child element but comma delimited
# This is preferred since file size is greatly reduced
#
# PARAMS
# input_filename = Input TSV file that will be given to read_file() (read documentation over there)
# output_filename = the output XML file
# is_prettify = Boolean representing whether the XML file should be prettified
# gene_count = limit of how many genes should be included in the XML. -1 is all genes in TSV
def create_xml_simple(input_filename, output_filename, is_prettify=True, gene_count=-1 ):
    payload = read_file(input_filename)
    data = payload[0]
    root = Element("genome")
    root.set('version', '1.0.0')
    root.append(Comment("Generated XML"))
    outer_count = 0;
    for entry in data:
        if gene_count != -1 and outer_count > gene_count:
            break
        if len(entry.rnaSeq) == 104:
            gene_child = SubElement(root, "gene")
            gene_db_id_child = SubElement(gene_child, "gene-db-id")
            gene_sec_id_child = SubElement(gene_child, "gene-sec-id")
            gene_name_child = SubElement(gene_child, "gene-name")
            gene_db_id_child.text = entry.dbIdentifier
            gene_sec_id_child.text = entry.secondaryIdentifier
            gene_name_child.text = entry.geneName
            rna_exp_child = SubElement(gene_child, "rna-exp")
            rna_exp_text = ""
            for rna_exp_stage in entry.rnaSeq:
                rna_exp_text += str(rna_exp_stage) + ","
            rna_exp_child.text = rna_exp_text[:-1]
            outer_count += 1
    file = open(output_filename, 'w')
    if is_prettify:
        file.write(prettify(root))
    else:
        file.write(str(tostring(root)))
    file.close()
                
from xml.etree import ElementTree
from xml.dom import minidom

#
# Input = Element Tree
# Intermediate process = Converts element tree to string, then uses mindom to prettyify string which becomes the output
# output = String prettified
#
def prettify(elem):
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


if __name__ == '__main__':
    # The input file contains all RNA expression data for all genes available in D. m
    create_xml_simple("../res/results_w_more_ids.tsv", "../res/genedata_w_ids.xml")
    #read_file("../res/results_w_more_ids.tsv")
        
             