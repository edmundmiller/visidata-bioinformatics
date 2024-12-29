'''Loaders for Bioinformatics file formats'''

__author__ = 'Edmund Miller <git@edmundmiller.dev>'
__version__ = '0.1'

vd.option('disp_hello', 'Hello world!', 'string to display for hello-world command')

BaseSheet.addCommand('0', 'hello-world', 'status(options.disp_hello)')
