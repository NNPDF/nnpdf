import os

import click

# https://github.com/pallets/click/blob/main/examples/complex/complex/cli.py
class ComplexCLI(click.MultiCommand):
    def list_commands(self, ctx):
        rv = []
        for filename in os.listdir(os.path.join(os.path.dirname(__file__), "commands")):
            if filename.endswith(".py") and filename.startswith("cmd_"):
                rv.append(filename[4:-3])
        rv.sort()
        return rv

    def get_command(self, ctx, name):
        try:
            mod = __import__(f"validphys.scripts.commands.cmd_{name}", None, None, ["cli"])
        except ImportError:
            return
        return mod.cli

@click.command(cls=ComplexCLI)
def main():
    """Here is some documentation
    """
    pass
