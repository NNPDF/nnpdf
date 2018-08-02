import sys
import os
import logging

#TODO: Look into making these lazy imports
import prompt_toolkit
from prompt_toolkit.contrib.completers import WordCompleter

from reportengine.colors import t

from validphys.app import App
from validphys import comparefittemplates

log = logging.getLogger(__name__)


class CompareFitApp(App):
    def add_positional_arguments(self, parser):
        parser.add_argument(
            'base_fit',
            default=None,
            nargs='?',
            help="The fit to produce the report for.",
        )
        parser.add_argument(
            'reference_fit',
            default=None,
            nargs='?',
            help="The fit to compare with")
        #These are not really positional, but it here they show up before others.
        parser.add_argument(
            '--title', help="The title that will be indexed with the report.")
        parser.add_argument('--author', help="The author of the report.")
        parser.add_argument(
            '--keywords', nargs='+', help="keywords to index the report with.")
        parser.add_argument(
            '-i',
            '--interactive',
            help="Ask interactively for the missing data",
            action='store_true')

    def try_complete_args(self):
        args = self.args
        argnames = ('base_fit', 'reference_fit', 'title', 'author', 'keywords')
        bad = [argname for argname in argnames if not args[argname]]
        if bad and not args['interactive']:
            sys.exit(f"The following arguments are required: {bad}")
        for arg in bad:
            self.args[arg] = getattr(self, f'interactive_{arg}')()
        texts = '\n'.join(
            f'    {argname.replace("_", " ").capitalize()}: {args[argname]}'
            for argname in argnames)
        log.info(f"Starting NNPDF fit comparison:\n{texts}")

    def interactive_base_fit(self):
        l = self.environment.loader
        completer = WordCompleter(l.available_fits)
        return prompt_toolkit.prompt("Enter base fit: ", completer=completer)

    def interactive_reference_fit(self):
        l = self.environment.loader
        completer = WordCompleter(l.available_fits)
        return prompt_toolkit.prompt(
            "Enter reference fit: ", completer=completer)

    def interactive_title(self):
        #TODO Use the colors in prompt_toolkit 2+ instead of this
        default = (f"Comparison between {self.args['base_fit']} "
                   f"and {self.args['reference_fit']} ")
        print(f"Enter report title [default:\n{t.dim(default)}]:")
        #Do not use the default keyword because it is a pain to delete
        res = prompt_toolkit.prompt("")
        if not res:
            return default
        return res

    def interactive_author(self):
        default = ""
        try:
            import pwd
        except ImportError:
            pass
        else:
            default = pwd.getpwuid(os.getuid())[4]
        return prompt_toolkit.prompt("Enter author name: ", default=default)

    def interactive_keywords(self):
        kwinp = prompt_toolkit.prompt("Enter keywords: ")
        return [k.strip() for k in kwinp.split(',') if k]

    def get_commandline_arguments(self, cmdline=None):
        args = super().get_commandline_arguments(cmdline)
        # This is needed because the environment wants to know how to resolve
        # the relative paths to find the templates. Best to have the template
        # look as much as possible as a runcard passed from the command line
        args['config_yml'] = comparefittemplates.template_path
        return args

    def get_config(self):
        self.try_complete_args()
        return super().get_config()


def main():
    a = CompareFitApp()
    a.main()


if __name__ == '__main__':
    main()
