import logging
import os
import pwd

from reportengine.compat import yaml
from validphys import deltachi2templates
from validphys.app import App

log = logging.getLogger(__name__)


class HyperoptPlotApp(App):
    def add_positional_arguments(self, parser):
        """Wrapper around argumentparser"""
        parser.add_argument("fit", help="Name of the fit")
        parser.add_argument("hessian_pdfs", help="Name of the set of Hessian pdfs")
        parser.add_argument("--Q", help="Energy Scale in GeV", type=float, default=1.7)
        # Report meta data
        parser.add_argument(
            "--author",
            help="Add custom author name to the report's meta data",
            type=str,
            default=pwd.getpwuid(os.getuid())[4].replace(",", ""),
        )
        parser.add_argument("--title", help="Add custom title to the report's meta data", type=str)
        parser.add_argument(
            "--keywords",
            help="Add keywords to the report's meta data. The keywords must be provided as a list",
            type=list,
            default=[],
        )

    def complete_mapping(self):
        args = self.args

        fit = args["fit"]
        hessian_pdfs = args["hessian_pdfs"]
        if args["title"] == None:
            args["title"] = f"$\Delta \chi^2$ report for {fit}"

        autosettings = {}
        autosettings["meta"] = {
            "title": args["title"],
            "author": args["author"],
            "keywords": args["keywords"],
        }
        autosettings["Q"] = args["Q"]
        autosettings["fit"] = fit
        autosettings["pdfs"] = [fit]
        autosettings["hessianinfo"] = {
            "fit": fit,
            "pdf": hessian_pdfs,
            "theory": {"from_": "fit"},
            "theoryid": {"from_": "theory"},
            "use_cuts": "fromfit",
            "dataset_inputs": {"from_": "fit"},
            "normalize_to": fit,
        }

        autosettings["decomposition"] = {"normalize_to": hessian_pdfs, "pdf": hessian_pdfs}
        autosettings["MC_Hessian_compare"] = {"pdfs": [hessian_pdfs, fit], "normalize_to": fit}

        return autosettings

    def get_config(self):
        complete_mapping = self.complete_mapping()
        runcard = deltachi2templates.template_path
        # No error handling here because this is our internal file
        with open(runcard) as f:
            # TODO: Ideally this would load round trip but needs
            # to be fixed in reportengine.
            c = yaml.safe_load(f)
        c.update(complete_mapping)
        return self.config_class(c, environment=self.environment)


def main():
    app = HyperoptPlotApp()
    app.main()


if __name__ == "__main__":
    main()
