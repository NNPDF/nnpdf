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
stub_path.mkdir(exist_ok=True)

this_path = pathlib.Path(__file__).parent

def gen_stubs():
    pattern = r'\[\[module\]\]'
    with (this_path/'swig.template').open() as f:
        template = f.read()
    for module in modules:
        stub = re.sub(pattern,module,template)
        p = stub_path / (module + '.i') 
        with p.open('w') as f:
            f.write(stub)

if __name__ == '__main__':
    gen_stubs()
