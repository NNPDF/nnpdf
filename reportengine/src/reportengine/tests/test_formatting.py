# -*- coding: utf-8 -*-
"""
Created on Mon May 16 14:49:52 2016

@author: Zahari Kassabov
"""
import unittest

from reportengine import namespaces, formattingtools


class Dataset:
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return "Nice: %s" % self.name

l = namespaces.NSList([Dataset("Hello"), Dataset("World")], nskey='dataset', )



class TestFormatting(unittest.TestCase):

    def setUp(self):
        self.ns = namespaces.ChainMap({'a': {'say': "Hi!"}, 'dataset':l})

    def test_formatting(self):
        name = formattingtools.get_nice_name(self.ns, ('a',('dataset', 0),))
        self.assertEqual(name, 'a_NiceHello')


if __name__ == '__main__':
    unittest.main()
