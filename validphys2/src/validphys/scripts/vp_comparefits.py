import sys
import os
import logging

#TODO: Look into making these lazy imports
import prompt_toolkit
from prompt_toolkit.completion import WordCompleter

from reportengine.compat import yaml
from reportengine.colors import t

from validphys.app import App
from validphys.loader import RemoteLoader
from validphys import comparefittemplates, compareclosuretemplates
from validphys.promptutils import confirm, KeywordsWithCache

log = logging.getLogger(__name__)

CURRENT_FIT_LABEL_DEFAULT = "Current Fit"
REFERENCE_FIT_LABEL_DEFAULT = "Reference Fit"



class CompareFitApp(App):
    def add_positional_arguments(self, parser):
        parser.add_argument(
            'current_fit',
            default=None,
            nargs='?',
            help="The fit to produce the report for.",
        )
        parser.add_argument(
            'reference_fit',
            default=None,
            nargs='?',
            help="The fit to compare with.")
        # Group together mandatory arguments that are not positional
        mandatory = parser.add_argument_group("mandatory", "Mandatory command line arguments")
        mandatory.add_argument(
            '--title', help="The title that will be indexed with the report.")
        mandatory.add_argument('--author', help="The author of the report.")
        mandatory.add_argument(
            '--keywords', nargs='+', help="keywords to index the report with.")

        parser.add_argument(
            '--thcovmat_if_present',
            action='store_true',
            help="Use theory cov mat for calculating statistical estimators if available.")
        parser.add_argument(
            '--no-thcovmat_if_present',
            action='store_true',
            help="DEPRECATED: does nothing")
        parser.add_argument(
            '--current_fit_label',
            nargs='?',
            default=CURRENT_FIT_LABEL_DEFAULT,
            help="The label for the fit that the report is being produced for.",
        )
        parser.add_argument(
            '--reference_fit_label',
            nargs='?',
            default=REFERENCE_FIT_LABEL_DEFAULT,
            help="The label for the fit that is being compared to.")
        parser.add_argument(
            '--dataset_join',
            nargs='?',
            default="intersection",
            choices=['intersection', 'union', 'current', 'reference'],
            help="Which datasets to use for the computation of chi2 for the data-theory comparison.")
        parser.add_argument(
            '-i',
            '--interactive',
            help="Ask interactively for the missing data",
            action='store_true')
        parser.add_argument(
            '-c',
            '--closure',
            help="Use the closure comparison template.",
            action='store_true')

        parser.set_defaults()

    def try_complete_args(self):
        args = self.args
        argnames = (
            'current_fit', 'reference_fit', 'title', 'author', 'keywords')
        optionalnames = (
            'current_fit_label', 'reference_fit_label', 'dataset_join')
        boolnames = (
            'thcovmat_if_present',)
        badargs = [argname for argname in argnames if not args[argname]]
        badbools = [bname for bname in boolnames if args[bname] is None]
        bad = badargs + badbools
        if bad and not args['interactive']:
            sys.exit(f"The following arguments are required: {bad}")
        try:
            for arg in bad:
                self.args[arg] = getattr(self, f'interactive_{arg}')()
            if args['interactive']:
                for arg in optionalnames:
                    self.args[arg] = getattr(self, f'interactive_{arg}')()
        except EOFError:
            raise KeyboardInterrupt()
        texts = '\n'.join(
            f'    {argname.replace("_", " ").capitalize()}: {args[argname]}'
            for argname in [*argnames, *optionalnames, *boolnames])
        log.info(f"Starting NNPDF fit comparison:\n{texts}")

    def interactive_current_fit(self):
        l = self.environment.loader
        completer = WordCompleter(l.available_fits)
        return prompt_toolkit.prompt("Enter current fit: ", completer=completer)

    def interactive_current_fit_label(self):
        #TODO Use the colors in prompt_toolkit 2+ instead of this
        default = CURRENT_FIT_LABEL_DEFAULT
        print(f"Enter label for current fit [default:\n{t.dim(default)}]:")
        #Do not use the default keyword because it is a pain to delete
        res = prompt_toolkit.prompt("")
        if not res:
            return default
        return res

    def interactive_reference_fit(self):
        l = self.environment.loader
        completer = WordCompleter(l.available_fits)
        return prompt_toolkit.prompt(
            "Enter reference fit: ", completer=completer)

    def interactive_reference_fit_label(self):
        #TODO Use the colors in prompt_toolkit 2+ instead of this
        default = REFERENCE_FIT_LABEL_DEFAULT
        print(f"Enter label for reference fit [default:\n{t.dim(default)}]:")
        #Do not use the default keyword because it is a pain to delete
        res = prompt_toolkit.prompt("")
        if not res:
            return default
        return res

    def interactive_title(self):
        #TODO Use the colors in prompt_toolkit 2+ instead of this
        default = (f"Comparison between {self.args['current_fit']} "
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
        if isinstance(self.environment.loader, RemoteLoader):
            completer = WordCompleter(words=KeywordsWithCache(self.environment.loader))
        else:
            completer = None
        kwinp = prompt_toolkit.prompt(
            "Enter keywords: ",
            completer=completer,
            complete_in_thread=True,
        )
        return [k.strip() for k in kwinp.split(',') if k]

    def interactive_thcovmat_if_present(self):
        """Interactively fill in the `use_thcovmat_if_present` runcard flag. Which is True by default
        """
        message = ("Do you want to use the theory covariance matrix, if available,\n"
                   "to calculate the statistical estimators? ")
        return confirm(message, default=True)

    def interactive_dataset_join(self):
        """Interactively fill in the `dataset_join` runcard flag. Which is "intersection" by default
        """
        message = ("For the data theory comparsion, do you wish to use the datasets of: \n"
                   "current, reference, intersection, or union")
        return confirm(message, default="intersection")

    def get_commandline_arguments(self, cmdline=None):
        args = super().get_commandline_arguments(cmdline)
        # This is needed because the environment wants to know how to resolve
        # the relative paths to find the templates. Best to have the template
        # look as much as possible as a runcard passed from the command line
        if not args['closure']:
            args['config_yml'] = comparefittemplates.template_path
        else:
            #This doesn't print anything
            log.info(f"using closure test template.")
            args['config_yml'] = compareclosuretemplates.template_path
        return args

    def complete_mapping(self):
        args = self.args
        autosettings = {}
        autosettings['meta'] = {
            'title': args['title'],
            'author': args['author'],
            'keywords': args['keywords']
        }
        currentmap = {'id': args['current_fit'], 'label': args['current_fit_label']}
        autosettings['current'] = {
            'fit': currentmap,
            'pdf': currentmap,
            'theory': {
                'from_': 'fit'
            },
            'theoryid': {
                'from_': 'theory'
            },
            'speclabel': args['current_fit_label']
        }
        refmap = {'id': args['reference_fit'], 'label': args['reference_fit_label']}
        autosettings['reference'] = {
            'fit': refmap,
            'pdf': refmap,
            'theory': {
                'from_': 'fit'
            },
            'theoryid': {
                'from_': 'theory'
            },
            'speclabel': args['reference_fit_label']
        }
        autosettings['use_thcovmat_if_present'] = args['thcovmat_if_present']
        autosettings['data_theory_dataset'] = {'dataset_join': args['dataset_join'], "use_cuts": "internal"}
        return autosettings


    def get_config(self):
        self.try_complete_args()
        #No error handling here because this is our internal file
        with open(self.args['config_yml']) as f:
            #TODO: Ideally this would load round trip but needs
            #to be fixed in reportengine.
            c = yaml.safe_load(f)
        c.update(self.complete_mapping())
        return self.config_class(c, environment=self.environment)


def main():
    a = CompareFitApp()
    a.main()


if __name__ == '__main__':
    main()
