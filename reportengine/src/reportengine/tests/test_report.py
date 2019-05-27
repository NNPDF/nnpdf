#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 10:19:59 2017

@author: zah
"""
import pytest

from reportengine.resourcebuilder import ResourceBuilder
from reportengine.report import Config
from reportengine.environment import Environment
from reportengine import namespaces
import reportengine.report
from reportengine.templateparser import CustomParsingError, get_targets_and_replace

expected_parsed = (
"""This is a test with title My report
===================================

ñ


### Another title, inside l: First
Yet world stays the same:
ñ

But can get repeated:
ñ
ñ
Also multiple times on the same line: ñ, ñ!

Can nest even more:

####Nested
See?

####Nested
See?


This was the previous title: First


### Another title, inside l: Second
Yet world stays the same:
ñ

But can get repeated:
ñ
ñ
Also multiple times on the same line: ñ, ñ!

Can nest even more:

####Nested
See?

####Nested
See?


This was the previous title: Second



Can nest providers:

Processed First

Processed Second


And this the original title My report.

Done."""
                   )

expected_second = (
"""Second: My report"""
)

class Providers:
    def processed(self, title):
        return "Processed " + title

def test_processing(tmpdir):
    inp = {'template': 'test_new.md', 'config_rel_path':'.',
           'report_title':{
                'title': "My report",
            },
            'othertemplate': {'template': 'test_second.md'},
            'title': "Bad title",

           'world': "ñ",
           'l': [{'title': "First"}, {'title': "Second"}],
           'nested': {'title': "Nested"},

           }

    spec = ('report_title',)
    otherspec = ('othertemplate', 'report_title')
    rb = ResourceBuilder(Config(inp), [reportengine.report, Providers()],
                         fuzzytargets=[('template_text', spec, (), ()),
                                       ('template_text', otherspec, (), ()),

                                       ],
                         environment=Environment(output=str(tmpdir)))
    rb.resolve_fuzzytargets()
    rb.execute_sequential()


    res = namespaces.resolve(rb.rootns, spec)['template_text']
    assert res== expected_parsed

    otherres = namespaces.resolve(rb.rootns, otherspec)['template_text']
    assert otherres== expected_second

def test_bad_matches():
    from io import StringIO
    with pytest.raises(CustomParsingError):
        list(get_targets_and_replace(StringIO('{@a::b@}')))

    with pytest.raises(CustomParsingError):
        list(get_targets_and_replace(StringIO('{@ab^^xx@}')))

    with pytest.raises(CustomParsingError):
        list(get_targets_and_replace(StringIO('{@aa bb bad@}')))
