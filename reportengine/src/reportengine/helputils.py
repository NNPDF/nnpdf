# -*- coding: utf-8 -*-
"""
Utils for printing documentation and formatting it properly.

Created on Fri Jul  1 11:40:46 2016

@author: Zahari Kassabov
"""
import textwrap
import itertools
import inspect
import re
from collections import OrderedDict
from io import StringIO


from reportengine.utils import get_providers
from reportengine.colors import t


_highlight_map = {}
_highlight_cycle = itertools.cycle([t.green, t.cyan ,t.yellow, t.magenta, t.blue])
def get_highlight_color(highlight):
    if highlight not in _highlight_map:
        _highlight_map[highlight] = next(_highlight_cycle)
    return _highlight_map[highlight]

def sane_wrap(txt, width=70, initial_indent='    ', subsequent_indent='  '):
    """Wrap lines so that the occupy at most ``width`` characters
    (except for long words), in such a way that a single continuation
    newline is ignored, but two newlines (i.e. a paragraph break)
    are respected, as well as indented text."""
    #Remove leading and trailing newlines
    txt = txt.strip('\n')
    #Remove continuation whitespaces
    txt = re.sub(r'(.)\n([^\s])',r'\1 \2',txt)
    def _indent():
        yield initial_indent
        while True:
            yield subsequent_indent
    indent = _indent()

    def wraplines(txt):
        for line in txt.splitlines(keepends=True):
            ni = next(indent)
            line = f"{ni}{line}"
            if len(line) >= width:
                #Break by the last space before the with limit
                break_space = line.rfind(' ', len(ni), width)
                if break_space == -1:
                    #If not possible, break by the first space afer the limit
                    break_space = line.find(' ', width)
                    if break_space == -1:
                        yield line
                        continue
                yield line[:break_space] + '\n'
                #+1 removes the marching space
                yield from wraplines(line[break_space+1:])
            else:
                yield line
    return list(wraplines(txt))

def sane_fill(txt, *args, **kwargs):
    return ''.join(sane_wrap(txt, *args, **kwargs))

def sane_dedent(txt):
    """Leave the first line alone and dedent the rest"""
    first_line = txt.find('\n')
    if first_line == -1:
        return txt
    else:
        return txt[:first_line+1] + textwrap.dedent(txt[first_line+1:])

def wrap_lines_with_color_header(color_header, white_header, rest,
                                 *,initial_indent='    ',  **kwargs):

    #This very complicated ligic is so that textwrap does not take into
    #account the length of the color encoding sequences when formatting.
    txt = "%s%s" % (white_header, sane_dedent(rest))

    #The default of removing double newlines is kind of nonsensical, so we
    #break lines manually


    txt = sane_fill(txt, initial_indent=initial_indent, **kwargs)
    txt = (initial_indent + color_header
           + txt[len(white_header) + len(initial_indent):])
    return txt

def print_signature(function):
    sig = inspect.signature(function)
    header = function.__name__
    color_header = t.bold(header)

    return wrap_lines_with_color_header(color_header, header, str(sig),
                                        initial_indent='')


def get_parser_type(f, sig_index=1):
    """Get a string corresponding to the valid type of a parser function.
    sig_index should be zero for instances and 1 for classes."""
    sig = inspect.signature(f)

    try:
        param_tp = list(sig.parameters.values())[sig_index].annotation
    except IndexError:
        return ''
    return get_annotation_string(param_tp)


def get_annotation_string(param_tp):
    if param_tp is inspect.Signature.empty:
        return ''

    if isinstance(param_tp, tuple):
        s = " or ".join(str(getattr(k,'__name__', k)) for k in param_tp)
    else:
        s = getattr(param_tp, '__name__' ,param_tp)
    return '(%s)' % s


def format_config_line(val, function, sig_index=1):
    #Get the docs
    doc = function.__doc__
    if doc is None:
        doc = ''

    #Get the recognized type
    tp = get_parser_type(function, sig_index=sig_index)



    color_header = "%s%s: " %(t.bold(val), tp)
    white_header = "%s%s: " %(val, tp)
    return wrap_lines_with_color_header(color_header, white_header, doc)


def format_environment(environ_class):
    header = t.bold_underline("Environment")+'\n\nThe following keys are injected from the environemnt:\n'
    lines = []
    for name, description in environ_class.ns_dump_description().items():
        expl = ": " + description
        line = wrap_lines_with_color_header(t.bold(name), name, expl)
        lines.append(line)
    return header+'\n\n'.join(lines)


def format_config(config_class):

    res = StringIO()

    all_parsers = config_class.get_all_parse_functions()
    header = (f"{t.bold_underline('Configuration parser')}"
              )

    config_header = f"{t.bold('Input')}\nThe following keys of the config file have a special meaning:\n"

    lines = []

    for val, function in all_parsers.items():

        lines.append(format_config_line(val, function))

    all_producer = config_class.get_all_produce_functions()
    produce_header = (f"{t.bold('Production rules')}\nThe following elements "
        "can be obtained as a function "
        "of the configuration keys:\n")

    produce_lines = []
    for val, function in all_producer.items():
        produce_lines.append(format_config_line(val, function))


    res.write(header)
    if lines:
        res.write(f'\n\n{config_header}')
        res.write('\n\n'.join(lines))
    if produce_lines:
        res.write(f'\n\n{produce_header}')
        res.write('\n\n'.join(produce_lines))
    return res.getvalue()

def print_providertree(providertree, environ_class=None):
    it = iter(providertree)
    result = StringIO()

    seen = set()
    res_providers = {}
    unknown_lines = OrderedDict()
    config_lines = OrderedDict()

    function = next(it)

    result.write(t.bold_underline(function.__name__))
    result.write('\n\n')
    result.write("Defined in: %s\n\n" % t.blue_underline(function.__module__))


    if hasattr(function, 'highlight'):
        hl = function.highlight
        result.write("Generates: %s\n\n" % get_highlight_color(hl)(hl))

    result.write(print_signature(function))
    result.write('\n\n')

    doc = function.__doc__
    if doc is not None:
        result.write(sane_fill(sane_dedent(doc),
                               initial_indent='', subsequent_indent=''))
        result.write('\n\n')


    def walk(it, provider=''):
        for spec in it:
            tp, name, value = spec
            #If we require it for the current provider '' drop the past one.
            if name not in res_providers or res_providers[name] != '':
                res_providers[name] = provider

            if name in seen:
                continue
            seen.add(name)
            if tp in ('config', 'produce'):
                if tp=='config':
                    config_lines[name] = format_config_line(name, value[0], sig_index=0)
                walk(value[1:], name)


            elif tp == 'provider':
                #We only care about the first level of nested providers.
                if provider=='':
                    walk(value[1:], name)
                else:
                    walk(value[1:], provider)

            elif tp == 'unknown':
                val_tp = get_annotation_string(value.annotation)

                if environ_class and name in environ_class.ns_dump_description():
                    default = t.bright_blue(" [Set based on the environment]")
                elif value.default is not value.empty:
                    default = ' = {}'.format(value.default)
                else:
                    default = ''
                line = "  {}{}{}".format(t.bold(name), val_tp, default)
                unknown_lines[name] = line
            else:
                raise ValueError("Unknown walk spec")

    walk(it)
    for reosurce,provider in res_providers.items():
        if provider:
            s = t.blue('[Used by %s' % t.underline(provider)) + t.blue(']')
            if reosurce in config_lines:
                config_lines[reosurce] += '\n  ' + s
            if reosurce in unknown_lines:
                unknown_lines[reosurce] += ' ' + s


    if config_lines:
        result.write("The following resources are read from the configuration:\n\n")
        #Sort so that the direct dependencies come first. Because the sort is
        #stable, the total ordering looks good this way.
        result.write('\n\n'.join(config_lines[k] for k in
                          sorted(config_lines, key=lambda x:res_providers[x]!='')))

    if unknown_lines:
        result.write(
        "\n\nThe following additionl arguments can "
        "be used to control the\nbehaviour. "
        "They are set by default to sensible values:\n\n"""
        )
        result.write('\n'.join(unknown_lines[k] for k in
                          sorted(unknown_lines, key=lambda x:res_providers[x]!='')))
    return result.getvalue()


def format_providermodule(module):
    moddoc = module.__doc__
    if moddoc is None:
        moddoc = ''
    moddoc = sane_fill(moddoc, initial_indent='', subsequent_indent='')

    functions = get_providers(module)
    lines = []

    name = t.bold_underline(module.__name__)

    for val, function in functions.items():

        #Get the docs
        doc = function.__doc__
        if doc is None:
            doc = ''

        if hasattr(function, 'highlight'):
            highlight = '(%s)' % function.highlight
            color_highlight = get_highlight_color(function.highlight)(highlight)
        else:
            color_highlight = highlight = ''


        color_header = "%s%s: "%(t.bold(val), color_highlight)
        white_header = "%s%s:"%(val, highlight)
        lines.append(wrap_lines_with_color_header(color_header, white_header, doc))



    lines = '\n\n'.join(lines)



    s = ("{name}\n{moddoc}\n"
        "The following providers are defined in this module:\n\n"
        "{lines}".format(name=name, moddoc = moddoc, lines=lines))
    return s
