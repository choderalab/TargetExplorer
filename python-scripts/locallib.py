from lxml import etree
from lxml.builder import E

class PDBCnstrct_Webpage(object):

    def __init__(self):
        self.html_tree = self.init_html_tree()
        self.head_node = self.html_tree.find('head')
        self.div_container = self.html_tree.find('body/div[@id="container"]')
        #self.main_content_node = self.html_tree.find('body/div/div/section[@id="main_content"]')

    def init_html_tree(self):
        html_node = E.html()
        head_node = etree.SubElement(html_node, 'head')
        body_node = etree.SubElement(html_node, 'body')
        div_container = etree.SubElement(body_node, 'div')
        div_container.set('id', 'container')
        #div_inner = etree.SubElement(div_container, 'div')
        #div_inner.set('class', 'inner')

        return html_node

    def add_standard_header(self):
        meta_node1 = etree.SubElement(self.head_node, 'meta')
        meta_node1.set('charset', 'utf-8')
        meta_node2 = etree.SubElement(self.head_node, 'meta')
        meta_node2.set('http-equiv', 'X-UA-Compatible')
        meta_node2.set('content', 'chrome=1')

        link_node1 = etree.SubElement(self.head_node, 'link')
        link_node1.set('href', 'https://fonts.googleapis.com/css?family=Chivo:900')
        link_node1.set('rel', 'stylesheet')
        link_node1.set('type', 'text/css')
        link_node2 = etree.SubElement(self.head_node, 'link')
        link_node2.set('rel', 'stylesheet')
        link_node2.set('type', 'text/css')
        link_node2.set('href', 'stylesheets/stylesheet.css')
        link_node2.set('media', 'screen')
        link_node3 = etree.SubElement(self.head_node, 'link')
        link_node3.set('rel', 'stylesheet')
        link_node3.set('type', 'text/css')
        link_node3.set('href', 'stylesheets/pygment_trac.css')
        link_node3.set('media', 'screen')
        link_node4 = etree.SubElement(self.head_node, 'link')
        link_node4.set('rel', 'stylesheet')
        link_node4.set('type', 'text/css')
        link_node4.set('href', 'stylesheets/print.css')
        link_node4.set('media', 'print')

        title_node = etree.SubElement(self.head_node, 'title')
        title_node.text = 'TargetExplorerDB by choderalab'

        header_node = etree.SubElement(self.div_container, 'header')
        header_node.set('class', 'TEDBheader')
        h1_node = etree.SubElement(header_node, 'h1')
        h1_node.text = 'TargetExplorerDB'
        #h2_node = etree.SubElement(header_node, 'h2')
        #h2_node.text = ''

    def add_website_footer(self):
        pass

    def add_content(self, content_node):
        self.div_container.append(content_node)

    def write_html(self, ofilename):
        ofile = open(ofilename, 'w')
        ofile.write( '<!DOCTYPE html>\n' + etree.tostring(self.html_tree, pretty_print=True) )
        ofile.close()

def gen_main_content_node():
    main_content_node = E.section()
    main_content_node.set('id', 'main_content')
    return main_content_node

def gen_PDBCnstrct_section_node():
    PDBCnstrct_content_node = E.section()
    PDBCnstrct_content_node.set('class', 'PDBCnstrctMain')
    return PDBCnstrct_content_node

