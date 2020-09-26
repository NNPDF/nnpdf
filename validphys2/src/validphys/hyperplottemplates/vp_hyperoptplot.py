from validphys.app import App
from validphys import hyperplottemplates
from reportengine.compat import yaml
import pwd
import os


class HyperoptPlotApp(App):
    def add_positional_arguments(self, parser):
        """ Wrapper around argumentparser """
        # Hyperopt settings
        parser.add_argument(
            "hyperopt_folder",
            help="Folder of the hyperopt fit to generate the report for",
        )
        parser.add_argument(
            "-l",
            "--loss_target",
            help="Choice for the definition of target loss",
            choices=['average', 'best_worst', 'std'],
            default='average',
        )
        parser.add_argument(
            "-v",
            "--val_multiplier",
            help="Fraction to weight the validation loss with (test_multiplier = 1-val_multiplier)",
            type=float,
            default=0.0,
        )
        parser.add_argument(
            "-if",
            "--include_failures",
            help="Flag to include failed runs in  the plots",
            action="store_true",
        )
        parser.add_argument(
            "-t",
            "--threshold",
            help="Value of the loss function from which to consider a run to have failed",
            type=float,
            default=1e3,
        )
        parser.add_argument(
            "-f",
            "--filter",
            help="Add the filter key=value to the dataframe",
            nargs="+",
            default=(),
        )
        parser.add_argument(
            "-c",
            "--combine",
            help="If more than one replica folder is found, combine all trials",
            action="store_false",
        )
        # Autofiltering
        parser.add_argument(
            "--autofilter",
            help="Given a number of keys, perform an autofilter (removing combinations of elements with worse rewards",
            nargs="+",
        )
        # Debugging
        parser.add_argument("--debug", help="Print debug information", action="store_true")
        # Report meta data
        parser.add_argument(
            "--author",
            help="Add custom author name to the report's meta data",
            type=str,
            default=pwd.getpwuid(os.getuid())[4],
        )
        parser.add_argument(
            "--title",
            help="Add custom title to the report's meta data",
            type=str,
        )
        parser.add_argument(
            "--keywords",
            help="Add keywords to the report's meta data. The keywords must be provided as a list",
            type=list,
            default=[]
        )
        args = parser.parse_args()
        return args

    def complete_mapping(self):
        args = self.args

        while args["hyperopt_folder"][-1] == "/":
            args["hyperopt_folder"] = args["hyperopt_folder"][:-1]
        folder_path = args["hyperopt_folder"]
        index_slash = folder_path.rfind("/") + 1
        name_folder = folder_path[index_slash:]
        args["title"] = f"NNPDF hyperoptimization report for {name_folder}"

        autosettings = {}
        autosettings["meta"] = {
            "title": args["title"],
            "author": args["author"],
            "keywords": args["keywords"],
        }
        autosettings["commandline_args"] = {
            "hyperopt_folder": args["hyperopt_folder"],
            "val_multiplier": args["val_multiplier"],
            "include_failures": args["include_failures"],
            "threshold": args["threshold"],
            "filter": args["filter"],
            "combine": args["combine"],
            "autofilter": args["autofilter"],
            "debug": args["debug"],
            "loss_target": args["loss_target"]
        }

        return autosettings

    def get_config(self):
        runcard = hyperplottemplates.template_path
        # No error handling here because this is our internal file
        with open(runcard) as f:
            # TODO: Ideally this would load round trip but needs
            # to be fixed in reportengine.
            c = yaml.safe_load(f)
        c.update(self.complete_mapping())
        return self.config_class(c, environment=self.environment)


def main():
    app = HyperoptPlotApp()
    app.main()


if __name__ == "__main__":
    main()
