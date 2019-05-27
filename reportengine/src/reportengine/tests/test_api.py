# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 15:00:00 2019

@author: Zahari Kassabov, Michael Wilson
"""

import unittest

from reportengine import configparser, api

first_input = {
    'input_a': 'egg',
    'input_b': 'spam',
    'time': '8AM'
}

second_input = {
    'restaurant': 'La Patata'
}

bad_provider = 'invalidphys'

class Restaurant:
    def __init__(self, description):
        self.description = description
    def as_input(self):
        return {
            'description': self.description,
            'time': '10AM'}

class Config(configparser.Config):
    def parse_input_a(self, input_a):
        return f'boiled {input_a}'
    def parse_input_b(self, input_b):
        return input_b
    def parse_restaurant(self, restaurant):
        return Restaurant(restaurant)
    def produce_restaurant_time(self):
        return self.parse_from_('restaurant', 'time', write=False)[1]

class TestEnv:
    def __init__(self, **kwargs):
        self.kwargs = kwargs

class Providers:

    def processed_spam(self):
        return "processed spam"

    def breakfast(self, input_a, input_b='ham', time="11AM"):
        return f"breakfast of {input_a} and {input_b} at {time}."

    class Juice:
        def __init__(self):
            pass


class TestAPI(unittest.TestCase):
    def test_badprovider(self):
        with self.assertRaises(ImportError):
            api.API([bad_provider], Config, TestEnv)

    def test_API_parsesinput(self):
        a = api.API([Providers()], Config, TestEnv)
        self.assertIsInstance(a.restaurant(restaurant="whetherspoons"), Restaurant)

    def test_API_getattr(self):
        a = api.API([Providers()], Config, TestEnv)
        self.assertIsInstance(a.Juice(), Providers.Juice)

    def test_API_simple_input(self):
        a = api.API([Providers()], Config, TestEnv)
        assert a.breakfast(**first_input) == "breakfast of boiled egg and spam at 8AM."

    def test_API_production_rule(self):
        a = api.API([Providers()], Config, TestEnv)
        assert a.restaurant_time(**second_input) == "10AM"

if __name__ == '__main__':
    unittest.main()
