from genelist_creator import read_file
from xml.etree.ElementTree import Element, SubElement, Comment, tostring

def create_xml(filename):
    payload = read_file(filename)
    data = payload[0]
    root = Element("genome")
    root.set('version', '1.0.0')
    root.append(Comment("Generated XML"))
    outer_count = 0;
    for entry in data:
        #if outer_count > 1000:
         #   break
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
    file = open("../res/data.xml", 'w')
    file.write(prettify(root))
    file.close()
                
from xml.etree import ElementTree
from xml.dom import minidom

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")
        
             