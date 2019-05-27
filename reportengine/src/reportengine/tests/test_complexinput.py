# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 21:09:40 2016

@author: Zahari Kassabov
"""
import unittest
from collections import ChainMap

import pytest

from reportengine import namespaces, configparser, utils, resourcebuilder
from reportengine.resourcebuilder import FuzzyTarget
from reportengine import collect
from reportengine.checks import make_argcheck

inp = {'pdfsets': ['a', 'b'],
       'theories': [1,2],
       'datasets': ['d1', 'd2'],
       'use_cuts': False,
       'cuts': {'use_cuts':True},
       'nocuts': {'use_cuts':False},
       'fits': ['NLO', 'NNLO'],

       'description': {'from_': 'fit'},
       'specialization': {
         'pdfsets': [{'from_': 'fit'}],
       },

       'maps': [
          {'fit': 'A',
           'pdfsets': ['X', {'from_': 'fit'}],
          },
          {'fit': 'B',
           'pdfsets': ['X', {'from_': 'fit'}],
          },
          {'fit': 'C',
           'pdfsets': ['X', {'from_': 'fit'}],
          },
        ],
        'ptos': [
          {'fit': 'X1',},
          {'fit': 'X2',},
        ],

       'fromeverywhere':{
               'fit':'N3LO',
               'pdf': 'XLO',
               'pdfsets':
               [
                       'XX',
                       {'from_':None},
                       {'from_':'fit'},
               ],

       },
       'datasepcs':[
               {'speclabel': 'l1'},
               {'nothing': True},
        ],
       't0spec':[

               {'use_t0': True,
                'pdf': 'T0PDF'
                },
               {'use_t0':False},

        ],
       'autons': {'namespaces_': 'nocuts::pdfsets::theories::datasets'},
       'nspiece': {'namespaces_': 'nocuts::pdfsets'},
       'badautons': {'namespaces_': 'nspiece::theories::datasets'},
}

other_input = {
    "fits": [{"from_": None}],
    "pdf": "not_used",
    "fit": "not_used",
    "dataspecs": [
        {
            "pdf": {"from_": "fit"},
            "dataspecs": [
                {"fit": "first_used", "pdfsets": [{"from_": None}]},
                {"fit": "second_used", "pdfsets": [{"from_": None}]},
            ],
        }
    ],
}

dataspecs_input = {
    "th": {"from_": "fit"},
    "use_cuts": True,
    "theory": {"from_": "th"},
    "pdf": {"from_": "fit"},
    "datasets": {"from_": "fit"},
    "dataspecs": [
        {"fit": "NNPDF31_nlo_as_0118_1000"},
        {"fit": "NNPDF31_nnlo_as_0118_1000"},
    ],
}

class Fit:
    def __init__(self, description):
        self.description = description

    def as_input(self):
        return {
            'description': self.description,
            'pdf': self.description,
            'datasets': [self.description, self.description, 'COMMON'],
            'th': {
                'theory': f'{self.description}'
            }
        }

class Config(configparser.Config):
    @configparser.element_of('pdfsets')
    def parse_pdf(self, pdf):
        return 'PDF: ' + pdf

    @configparser.element_of('theories')
    def parse_theory(self, theory):
        return 'th ' + str(theory)


    @configparser.element_of('datasets')
    def parse_dataset(self, ds,  theory, use_cuts):
        return f'ds: {ds} (theory: {theory}, cuts: {use_cuts})'

    def parse_template(self, template, rel_path):
        return template

    def produce_template_text(self, template):
        return template

    def parse_use_t0(self, use:bool, pdf=None):
        return use

    def produce_t0(self, use_t0, pdf=None):
        if use_t0:
            assert pdf
            return pdf
        else:
            return None

    def parse_experiment_input(self, inp:str):
        return inp

    def produce_experiment(self, experiment_input):
        return "experiment: " + experiment_input

    def produce_implicit_exp(self):
        return {'experiment': 'experiment: IMPLICIT'}

    @configparser.element_of('fits')
    def parse_fit(self, description):
        return Fit(description)

    def produce_fitpdf(self):
        return {'pdf': self.parse_from_('fit', 'pdf', write=False)[1]}

    def produce_matched_datasets_from_dataspecs(self, dataspecs):
        all_names = []
        for spec in dataspecs:
            with self.set_context(ns=self._curr_ns.new_child(spec)):
                _, datasets = self.parse_from_(
                    None, 'datasets', write=False)
                names = {ds.split()[1]: ds for ds in datasets}
                all_names.append(names)
        used_set = set.intersection(*(set(d) for d in all_names))

        res = []
        for k in used_set:
            inres = {'dataset_name': k}
            inner_spec_list = inres['dataspecs'] = []
            for ispec, spec in enumerate(dataspecs):
                d = ChainMap({
                    'dataset': all_names[ispec][k],
                }, spec)
                inner_spec_list.append(d)
            res.append(inres)
        res.sort(key=lambda x: (x['dataset_name']))
        return res


@make_argcheck
def bad_check(pdf):
    return pdf

class Providers:
    def report(self, template_text):
        return template_text

    def plot_a_pdf(self, pdf):
        return "PLOT OF " + str(pdf)

    dataspecs_speclabel = collect('speclabel', ('datasepcs',),
                                  element_default='label')
    @bad_check
    def bad_plot(self, pdf):
        return self.plot_a_pdf(pdf)



class TestSpec(unittest.TestCase):
    def test_nsexpand(self):
        spec = ('pdfsets', 'theories', 'datasets')
        c = Config(inp)
        ns = utils.ChainMap()
        specs = c.process_fuzzyspec(spec, ns=ns)
        self.assertEqual(len(specs), 8)
        datasets = [
              'ds: d1 (theory: th 1, cuts: False)',
              'ds: d2 (theory: th 1, cuts: False)',
              'ds: d1 (theory: th 2, cuts: False)',
              'ds: d2 (theory: th 2, cuts: False)',
              'ds: d1 (theory: th 1, cuts: False)',
              'ds: d2 (theory: th 1, cuts: False)',
              'ds: d1 (theory: th 2, cuts: False)',
              'ds: d2 (theory: th 2, cuts: False)']
        for spec, ds in zip(specs, datasets):
            self.assertEqual(namespaces.resolve(ns, spec)['dataset'], ds)

    def test_nsspec(self):
        c = Config(inp)
        spec = ('pdfsets', 'theories', 'datasets')
        targets = [
                   FuzzyTarget('dataset', spec+(), (), ()),
                   FuzzyTarget('dataset', ('cuts',)+spec, (), ()),
                   FuzzyTarget('dataset', ('nocuts',)+spec, (), ()),
                   FuzzyTarget('dataset', ('autons',), (), ()),
                  ]
        builder = resourcebuilder.ResourceBuilder(c, (), targets)
        builder.resolve_fuzzytargets()
        ns = namespaces.resolve(builder.rootns,
                  ('cuts', ('pdfsets',0), ('theories', 0), ('datasets', 0),))

        assert(ns['use_cuts']==True)
        assert(ns['dataset']=="ds: d1 (theory: th 1, cuts: True)" )
        ns = namespaces.resolve(builder.rootns,
                  ('nocuts', ('pdfsets',0), ('theories', 0), ('datasets', 0),))

        assert(ns['use_cuts']==False)
        assert(ns['dataset']=="ds: d1 (theory: th 1, cuts: False)" )

        nsfromauto = namespaces.resolve(builder.rootns,
            (('autons', 0),
             'nocuts', ('pdfsets',0), ('theories', 0), ('datasets', 0),))
        assert(ns['use_cuts'] == nsfromauto['use_cuts'])
        assert(ns['dataset'] == nsfromauto['dataset'])

        ns = namespaces.resolve(builder.rootns,
                  (('pdfsets',0), ('theories', 0), ('datasets', 0),))

        assert(ns['use_cuts']==False)
        assert(ns['dataset']=="ds: d1 (theory: th 1, cuts: False)" )


    def test_composed_namespaces_(self):
        c = Config(inp)
        spec = ('pdfsets', 'theories', 'datasets')
        targets = [
                   FuzzyTarget('dataset', spec+(), (), ()),
                   FuzzyTarget('dataset', ('cuts',)+spec, (), ()),
                   FuzzyTarget('dataset', ('nocuts',)+spec, (), ()),
                   FuzzyTarget('dataset', ('badautons',), (), ()),
                  ]
        builder = resourcebuilder.ResourceBuilder(c, (), targets)
        builder.resolve_fuzzytargets()
        ns = namespaces.resolve(builder.rootns,
            ('nocuts', ('pdfsets',1), ('theories', 1), ('datasets', 1),))

        nsfromauto = namespaces.resolve(builder.rootns,
            (('badautons', 7),))

        assert(ns['use_cuts'] == nsfromauto['use_cuts'])
        assert(ns['dataset'] == nsfromauto['dataset'])

    def test_gets_dependencies(self):
        inp = {
         'template': 'template',
         'mapping': {'template': 'othertemplate'},
         #'rel_path': 'cochambres'
        }
        c = Config(inp)
        targets =  [
                    FuzzyTarget('report', (), (), ()),
                    FuzzyTarget('report', ('mapping',), (), ())
                   ]
        builder = resourcebuilder.ResourceBuilder(c, Providers(), targets)
        builder.rootns.update({'rel_path':'examples'})
        builder.resolve_fuzzytargets()
        builder.execute_sequential()
        assert namespaces.resolve(builder.rootns, ('mapping',))['report'] == 'othertemplate'

    def test_iter_from(self):
        c = Config(inp)
        targets = [FuzzyTarget('description', ('fits',), (), ())]
        builder = resourcebuilder.ResourceBuilder(c, Providers(), targets)
        builder.resolve_fuzzytargets()
        builder.execute_sequential()
        s1 = [('fits', 0)]
        s2 = [('fits', 1)]
        assert namespaces.resolve(builder.rootns, s1)['description'] == 'NLO'
        assert namespaces.resolve(builder.rootns, s2)['description'] == 'NNLO'
        assert 'description' not in builder.rootns

    def test_iter_from_2(self):
        c = Config(inp)
        targets = [FuzzyTarget('pdfsets', ('ptos', 'specialization'), (), ())]
        builder = resourcebuilder.ResourceBuilder(c, Providers(), targets)
        builder.resolve_fuzzytargets()
        builder.execute_sequential()
        s1 = [('ptos', 0), 'specialization']
        s2 = [('ptos', 1), 'specialization']
        assert namespaces.resolve(builder.rootns, s1)['pdfsets'] == ['PDF: X1']
        assert namespaces.resolve(builder.rootns, s2)['pdfsets'] == ['PDF: X2']

    def test_nested_from(self):
        c = Config(inp)
        targets = [FuzzyTarget('plot_a_pdf', ('maps','pdfsets'), (), ())]
        builder = resourcebuilder.ResourceBuilder(c, Providers(), targets)
        builder.resolve_fuzzytargets()
        builder.execute_sequential()
        s0 = (('maps', 0), )
        s1 = (('maps', 1), )
        s2 = (('maps', 2), )

        assert namespaces.resolve(builder.rootns, s0)['pdfsets'] == ['PDF: X', 'PDF: A']
        assert namespaces.resolve(builder.rootns, s1)['pdfsets'] == ['PDF: X', 'PDF: B']
        assert namespaces.resolve(builder.rootns, s2)['pdfsets'] == ['PDF: X', 'PDF: C']

    def test_from_none(self):
        c = Config(inp)
        s = ('fromeverywhere',)
        targets = [FuzzyTarget('plot_a_pdf', s, (), ())]
        builder = resourcebuilder.ResourceBuilder(c, Providers(), targets)
        builder.resolve_fuzzytargets()
        builder.execute_sequential()
        assert namespaces.resolve(builder.rootns, s)['pdfsets'] == ['PDF: XX', 'PDF: XLO', 'PDF: N3LO']

    def test_complex_produce(self):
        c = Config(inp)
        s = ('t0spec',)
        targets = [FuzzyTarget('t0', s, (), ())]
        builder = resourcebuilder.ResourceBuilder(c, Providers(), targets)
        builder.resolve_fuzzytargets()
        builder.execute_sequential()
        assert namespaces.resolve(builder.rootns, [('t0spec',0)])['t0'] == 'PDF: T0PDF'
        assert namespaces.resolve(builder.rootns, [('t0spec',1)])['t0'] is None

    def test_produce_priority(self):
        c = Config(inp)
        s = ['implicit_exp']
        targets = [FuzzyTarget('experiment', s, (), ())]
        builder = resourcebuilder.ResourceBuilder(c, Providers(), targets)
        builder.resolve_fuzzytargets()
        builder.execute_sequential()
        assert namespaces.resolve(builder.rootns, s)['experiment'] == 'experiment: IMPLICIT'

    def test_default_collect(self):
        c = Config(inp)
        targets = [FuzzyTarget('dataspecs_speclabel', (), (), ())]
        builder = resourcebuilder.ResourceBuilder(c, Providers(), targets)
        builder.resolve_fuzzytargets()
        assert builder.rootns['dataspecs_speclabel'] == ['l1', 'label']

    def test_bad_check(self):
        c = Config(inp)
        targets = [FuzzyTarget('bad_plot', ('t0spec',), (), ())]
        builder = resourcebuilder.ResourceBuilder(c, Providers(), targets)
        with pytest.raises(TypeError):
            builder.resolve_fuzzytargets()

    def test_from_forward(self):
        c = Config(other_input)
        targets = [
            FuzzyTarget('pdfsets', ('dataspecs', 'dataspecs', 'fitpdf', 'fits'),
                        (), ())
        ]
        builder = resourcebuilder.ResourceBuilder(c, Providers(), targets)
        builder.resolve_fuzzytargets()
        builder.execute_sequential()
        assert namespaces.resolve(builder.rootns,
                                  (('dataspecs', 0), ('dataspecs', 0), 'fitpdf',
                                   ('fits', 0)))['pdfsets'][0] == 'PDF: first_used'
        assert namespaces.resolve(
            builder.rootns, (('dataspecs', 0), ('dataspecs', 1), 'fitpdf',
                             ('fits', 0)))['pdfsets'][0] == 'PDF: second_used'

    def test_namespace_production(self):
        c = Config(dataspecs_input)
        targets = [FuzzyTarget('dataset', ('matched_datasets_from_dataspecs', 'dataspecs'), (), ())]
        builder = resourcebuilder.ResourceBuilder(c, Providers(), targets)
        builder.resolve_fuzzytargets()
        builder.execute_sequential()
        ds1 = namespaces.resolve(builder.rootns,
                (('matched_datasets_from_dataspecs', 0), ('dataspecs', 0),))['dataset']
        ds2 = namespaces.resolve(builder.rootns,
                (('matched_datasets_from_dataspecs', 0), ('dataspecs', 1),))['dataset']
        assert ds1 != ds2

if __name__ == '__main__':
    unittest.main()
