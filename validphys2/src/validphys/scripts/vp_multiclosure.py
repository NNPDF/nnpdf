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
from validphys import compareinconsistentclosuretemplates
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
            '--lambdasettings_path',
            help="The path to the yaml file containing the settings for the lambdavalues report.")
        
        parser.set_defaults()

    def try_complete_args(self):
        args = self.args
        argnames = (
            'current_fit', 'reference_fit', 'title', 'author', 'keywords')
        optionalnames = (
            'current_fit_label', 'reference_fit_label')
        badargs = [argname for argname in argnames if not args[argname]]
        bad = badargs
        if bad:
            sys.exit(f"The following arguments are required: {bad}")
        texts = '\n'.join(
            f'    {argname.replace("_", " ").capitalize()}: {args[argname]}'
            for argname in [*argnames, *optionalnames])
        log.info(f"Starting NNPDF fit comparison:\n{texts}")

    def get_commandline_arguments(self, cmdline=None):
        args = super().get_commandline_arguments(cmdline)
        # This is needed because the environment wants to know how to resolve
        # the relative paths to find the templates. Best to have the template
        # look as much as possible as a runcard passed from the command line
        args['config_yml'] = compareinconsistentclosuretemplates.template_path
        return args

    def complete_mapping(self):
        args = self.args
        autosettings = {}
        shortsettings = {}
        autosettings['meta'] = {
            'title': args['title'],
            'author': args['author'],
            'keywords': args['keywords']
        }
        currentmap = {'id': args['current_fit'], 'label': args['current_fit_label']}
        shortsettings['current'] = {
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
        shortsettings['reference'] = {
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
        autosettings["compare_settings"] = shortsettings
        return autosettings

    def complete_lambdavaluesmapping(self):
        args = self.args
        # opening conf file
        with open(args['lambdasettings_path']) as f:
            settings = yaml.safe_load(f)
        autosettings = {}
        list_lambdas = []
        for lambdas in settings:
            lambdasetting = {}
            for set in ["variancepdf", "t0pdfset", "explained_variance_ratio", "label", "fits", "fit"]:
                lambdasetting[set] = lambdas[set]
            list_lambdas.append(lambdasetting)
        autosettings["lambdavalues"] = list_lambdas
        return autosettings
            
    def get_config(self):
        self.try_complete_args()
        #No error handling here because this is our internal file
        with open(self.args['config_yml']) as f:
            #TODO: Ideally this would load round trip but needs
            #to be fixed in reportengine.
            c = yaml.safe_load(f)
        complete_mapping = self.complete_mapping()
        complete_mapping.update(self.complete_lambdavaluesmapping())
        c["meta"] = complete_mapping["meta"]
        for lambdas, lambdas_mapping in zip(c["lambdavalues"], complete_mapping["lambdavalues"]):
            for set in ["variancepdf", "t0pdfset", "explained_variance_ratio", "label", "fits", "fit"]:
                lambdas[set] = lambdas_mapping[set]
        for set in ["current", "reference"]:
            c["compare_settings"][set] = complete_mapping["compare_settings"][set]
        return self.config_class(c, environment=self.environment)


def main():
    a = CompareFitApp()
    a.main()


if __name__ == '__main__':
    main()
