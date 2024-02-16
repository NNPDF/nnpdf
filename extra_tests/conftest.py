"""
This enables the use of a command line argument to set regenerate, which is used to autmate the
regeneration of the regression data.
"""


def pytest_addoption(parser):
    parser.addoption("--regenerate", action="store", default=False)


def pytest_generate_tests(metafunc):
    option_value = metafunc.config.option.regenerate
    if 'regenerate' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("regenerate", [option_value])
