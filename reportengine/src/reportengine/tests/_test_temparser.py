# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 15:14:42 2015

@author: zah
"""

import unittest

import jinja2

from reportengine.templateparser import Environment, TemplateRecordError
from reportengine.resourcebuilder import ResourceBuilder

class Resources:

    def resource(self):
        return "Resource"

    def otherresource(self, param=None):
        return "<%s>" % param

class TestStubs(unittest.TestCase):
    def setUp(self):
        env = Environment(loader=jinja2.PackageLoader('reportengine'))
        self.env = env
        self.test_template = env.get_template("test.md")

    def test_parser(self):
        with self.env.fetch_mode():
            self.test_template.render()
            t = self.env.from_string(r"{{f(arg, arg, kw=1)}}")
            with self.assertRaises(TemplateRecordError):
                t.render()
        otherresource = lambda param : param
        resource = "Resource"
        title = "Hello"
        normal_render = self.test_template.render(otherresource=otherresource,
                                        resource=resource, title=title)
        self.assertTrue(normal_render.startswith("This is a test with title Hello"))

    def test_with_builder(self):
        namespace = {"title": "Title"}

        with self.env.fetch_mode():
            self.test_template.render()

        targets = self.env.targets


        builder = ResourceBuilder(Resources(), targets ,namespace)
        builder.build_graph()
        builder.execute_sequential()
        results = [namespace[k] for k in builder.target_keys]
        with self.env.subs_mode(results):
            result = self.test_template.render()
        should_be = 'This is a test with title Title\n===================================\n\nResource\n\n\n\n<0>\n\n\n\n<1>\n\n\n\n<2>\n\n\n\n<3>\n\n\n\n<4>\n\n\n\n<None>\n<3>\n'
        self.assertEqual(result, should_be)

        gen = self.env.render_with_targets(self.test_template)
        targets = gen.send(None)
        namespace = {"title": "Title"}
        builder = ResourceBuilder(Resources(), targets ,namespace)
        builder.build_graph()
        builder.execute_sequential()
        results = [namespace[k] for k in builder.target_keys]
        result = gen.send(results)
        self.assertEqual(result, should_be)


if __name__ == "__main__":
    unittest.main()
