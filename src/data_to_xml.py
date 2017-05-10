from genelist_creator import read_file_rnadata
from genelist_creator import read_file_geneid
from xml.etree.ElementTree import Element, SubElement, Comment, tostring

# DEPRECATED: The secondary column and gene name was removed from the corresponding xml. Those columns are now read in the gene_id xml marshaller
#   
#    
# Creates XML from data where each rna stage is a new child element
def rna_exp_marshall_xml(input_filename, output_filename, is_prettify=True, gene_count=-1 ):
    payload = read_file_rnadata(input_filename)
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
def rna_exp_marshall_xml_simple(input_filename, output_filename, is_prettify=True, gene_count=-1 ):
    payload = read_file_rnadata(input_filename)
    data = payload[0]
    root = Element("genome-rna-data")
    root.set('version', '1.0.0')
    root.append(Comment("Generated XML by @CodingBash"))
    outer_count = 0;
    for entry in data:
        if gene_count != -1 and outer_count > gene_count:
            break
        if len(entry.rnaSeq) == 104:
            gene_child = SubElement(root, "gene")
            gene_db_id_child = SubElement(gene_child, "gene-db-id")
            gene_db_id_child.text = entry.dbIdentifier
            rna_exp_child = SubElement(gene_child, "rna-exp")
            rna_exp_text = ""
            for rna_exp_stage in entry.rnaSeq:
                rna_exp_text += str(rna_exp_stage) + ","
            rna_exp_child.text = rna_exp_text[:-1] # This removes the comma at the end
            outer_count += 1
    file = open(output_filename, 'w')
    if is_prettify:
        file.write(prettify(root))
    else:
        file.write(str(tostring(root)))
    file.close()
                
# Creates XML from data. This will marshall the the gene information (gene ID, secondary ID, name, and synonyms) to an XML
# The synonyms will be csv instead of separate childs to reduce file size
# PARAMS
# input_filename = Input TSV file that will be given to read_file() (read documentation over there)
# output_filename = the output XML file
# is_prettify = Boolean representing whether the XML file should be prettified
# gene_count = limit of how many genes should be included in the XML. -1 is all genes in TSV
def gene_id_marshall_xml_simple(input_filename, output_filename, is_prettify=True, gene_count=-1 ):
    payload = read_file_geneid(input_filename)
    data = payload[0]
    root = Element("genome-gene-ids")
    root.set('version', '1.0.0')
    root.append(Comment("Generated XML by @CodingBash"))
    outer_count = 0;
    for entry in data:
        if gene_count != -1 and outer_count > gene_count:
            break
        else:
            gene_child = SubElement(root, "gene")
            gene_db_id_child = SubElement(gene_child, "gene-db-id")
            gene_sec_id_child = SubElement(gene_child, "gene-sec-id")
            gene_name_child = SubElement(gene_child, "gene-name")
            gene_db_id_child.text = entry.dbIdentifier
            gene_sec_id_child.text = entry.secondaryIdentifier
            gene_name_child.text = entry.geneName
            synonyms_child = SubElement(gene_child, "synonyms")
            synonyms_text = ""
            for synonym in entry.synonyms:
                synonyms_text += str(synonym) + "\t" # tab separated instead of csv in case of intermediate commas
            synonyms_child.text = synonyms_text[:-1] # This removes the comma at the end
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
    
    rna_exp_marshall_xml_simple("../res/tsv_files/rna_data_limited_ids.tsv", "../res/xml_files/genome-rna-data-1.xml")
    gene_id_marshall_xml_simple("../res/tsv_files/gene_identifiers.tsv", "../res/xml_files/genome-gene-id-data-1.xml")
    
    #read_file("../res/results_w_more_ids.tsv")
        
             