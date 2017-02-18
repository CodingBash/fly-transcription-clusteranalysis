from genelist_creator import read_file
from xml.etree.ElementTree import Element, SubElement, Comment, tostring

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
                
from xml.etree import ElementTree
from xml.dom import minidom

def prettify(elem):
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")
        
             