# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 17:30:38 2016

@author: Zahari Kassabov
"""
import difflib

class ErrorWithAlternatives(Exception):

    alternatives_header="Instead of '%s', did you mean one of the following?"
    def __init__(self, message, bad_item = None, alternatives = None, *,
                 display_alternatives='best'):
        super().__init__(message)
        self.bad_item = bad_item
        if alternatives:
            alternatives = list(alternatives)
        self.alternatives = alternatives
        self.display_alternatives = display_alternatives

    def alternatives_text(self):
        if (self.display_alternatives=='none' or not self.display_alternatives
            or not self.alternatives):
            return ''
        if self.display_alternatives == 'best':
            alternatives = difflib.get_close_matches(self.bad_item,
                                                     self.alternatives)
        elif self.display_alternatives == 'all':
            alternatives = self.alternatives
        else:
            raise ValueError("Unrecognized display_alternatives option. "
            "Must be one of: 'all', 'best' or 'none'.")
        if not alternatives:
            return ''
        head = (self.alternatives_header
                % (self.bad_item,))
        txts = [' - {}'.format(alt) for alt in alternatives]
        return '\n'.join((head, *txts))

    def __str__(self):
        return '%s\n%s'%(self.args[0], self.alternatives_text())

class AsInputError(Exception): pass