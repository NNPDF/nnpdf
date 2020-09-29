import logging
import os
import pathlib
import pwd

import numpy as np

from reportengine.compat import yaml

from validphys import deltachi2templates, lhaindex
from validphys.api import API
from validphys.app import App
from validphys.core import PDF
from validphys.lhio import new_pdf_from_indexes

log = logging.getLogger(__name__)


class HyperoptPlotApp(App):
    def add_positional_arguments(self, parser):
        """ Wrapper around argumentparser """
        parser.add_argument(
            "fit",
            help="Name of the fit",
        )
        parser.add_argument(
            "hessian_pdfs",
            help="Name of the set of Hessian pdfs",
        )
        parser.add_argument(
            "--Q",
            help="Energy Scale in GeV",
            type=float,
            default=1.7,
        )
        parser.add_argument(
            "-t0",
            "--t0pdfset",
            help="PDF set used to generate the t0 covmat",
            type=str,
            default="NNPDF31_nnlo_as_0118",
        )
        # Report meta data
        parser.add_argument(
            "--author",
            help="Add custom author name to the report's meta data",
            type=str,
            default=pwd.getpwuid(os.getuid())[4].replace(',',''),
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
            default=[],
        )
        args = parser.parse_args()
        return args

    def complete_mapping(self):
        args = self.args

        fit = args["fit"]
        hessian_pdfs = args["hessian_pdfs"]
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
            "experiments": {"from_": "fit"},
            "use_t0": True,
            "t0pdfset": args["t0pdfset"],
        }

        autosettings["decomposition"] = {
            "pdfs": [hessian_pdfs, f"{hessian_pdfs}_pos", f"{hessian_pdfs}_neg"],
            "normalize_to": hessian_pdfs,
        }
        autosettings["MC_Hessian_compare"] = {
            "pdfs": [hessian_pdfs, fit],
            "normalize_to": fit,
        }

        return autosettings

    def get_config(self):
        def decompose(inputs):

            log.info("Decomposing Hessian pdfs...")

            pdf_name = inputs["pdf"]
            pdf = PDF(name=pdf_name)

            delta_chi2 = API.delta_chi2_hessian(**inputs)
            ind_pos = np.asarray([i for i in range(len(delta_chi2)) if delta_chi2[i] >= 0])
            ind_neg = np.asarray([i for i in range(len(delta_chi2)) if i not in ind_pos])

            ind_pos, ind_neg = ind_pos + 1, ind_neg + 1
            lhapdfpath = pathlib.Path(lhaindex.get_lha_datapath())
            new_pdf_from_indexes(
                pdf=pdf,
                indexes=ind_pos,
                folder=lhapdfpath,
                set_name=pdf_name + "_pos",
                hessian=True,
            )
            new_pdf_from_indexes(
                pdf=pdf,
                indexes=ind_neg,
                folder=lhapdfpath,
                set_name=pdf_name + "_neg",
                hessian=True,
            )

            log.info("Completed decomposing Hessian pdfs")

        complete_mapping = self.complete_mapping()

        decompose(complete_mapping["hessianinfo"])

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
