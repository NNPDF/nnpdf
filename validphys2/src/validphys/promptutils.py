# -*- coding: utf-8 -*-
"""Module which extends the functionality of promp_toolkit for user inputs/interactivity"""

# NOTE: since functions here were ported from wiki_upload, the uploads are as late as possible to
# help with command line speed
from prompt_toolkit import HTML

def yes_no_str(default=None):
    """Return a yes or no string for the prompt, with the default
    highlighted"""
    if default is None:
        return f'[y/n]'
    elif default:
        return HTML('[<b>Y</b>/n]')
    else:
        return HTML('[y/<b>N</b>]')

def confirm(message, default=None):
    """
    This is like prompt_toolkit.shortcuts.confirm (implemented by
    create_confirm_session) except that it doesn't bind control+c to "No", but
    instead raises an exception.

    It also support defaults.
    """
    from prompt_toolkit.key_binding.key_bindings import KeyBindings
    from prompt_toolkit.keys import Keys
    from prompt_toolkit.formatted_text import merge_formatted_text
    from prompt_toolkit.shortcuts import PromptSession

    bindings = KeyBindings()

    @bindings.add('y')
    @bindings.add('Y')
    def yes(event):
        session.default_buffer.text = 'y'
        event.app.exit(result=True)

    @bindings.add('n')
    @bindings.add('N')
    def no(event):
        session.default_buffer.text = 'n'
        event.app.exit(result=False)

    @bindings.add(Keys.Any)
    def nothing(event):
        " Disallow inserting other text. "
        pass

    if default:
        bindings.add(Keys.Enter)(yes)
    elif default is not None:
        bindings.add(Keys.Enter)(no)
    else:
        bindings.add(Keys.Enter)(nothing)

    complete_message = merge_formatted_text([message, yes_no_str(default)])
    session = PromptSession(complete_message, key_bindings=bindings)
    return session.prompt()

# We need some sort of cache because prompt_toolkit calls the callable
# every time it tries to complete.
class KeywordsWithCache:
    def __init__(self, loader):
        self.loader = loader
        self.words = None

    def __call__(self):
        if self.words is None:
            try:
                self.words = self.loader.remote_keywords
            # Catch a broad exception here as we don't want the completion
            # to break the app
            except Exception as e:
                self.words = []
        return self.words
