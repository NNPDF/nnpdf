import pathlib
import re

modules =[
 'logger',
 'commondata',
 'experiments',
 'parametrisation',
 'common',
 'utils',
 'fkset',
 'fkgenerator',
 'dataset',
 'lhapdfset',
 'randomgenerator',
 'pdfset',
 'fastkernel',
 'nnmpi',
 'positivity',
 'thpredictions',
 'timer',
 'nnpdfdb'
 ]

stub_path = pathlib.Path('./stubs')

interface_path = pathlib.Path('./src')

def gen_stubs():
    pattern = r'\[\[module\]\]'
    with open('swig.template') as f:
        template = f.read()
    for module in modules:
        stub = re.sub(pattern,module,template)
        p = stub_path / (module + '.i') 
        with p.open('w') as f:
            f.write(stub)

if __name__ == '__main__':
    gen_stubs()
