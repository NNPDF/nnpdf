from validphys.app import App
from validphys import deltachitemplates
from reportengine.compat import yaml
import pwd
import os
from validphys.core import PDF
from validphys.api import API
import numpy as np
from validphys.lhio import new_pdf_from_indexes
import logging

log = logging.getLogger(__name__)

class HyperoptPlotApp(App):
    def add_positional_arguments(self, parser):
        """ Wrapper around argumentparser """
        parser.add_argument(
            "original_pdfs",
            help="Name of the set of original pdfs",
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
            default=[],
        )
        args = parser.parse_args()
        return args

    def complete_mapping(self):
        args = self.args

        original_pdfs = args["original_pdfs"]
        hessian_pdfs = args["hessian_pdfs"]
        args["title"] = f"$\Delta \chi^2$ report for {original_pdfs}"

        autosettings = {}
        autosettings["meta"] = {
            "title": args["title"],
            "author": args["author"],
            "keywords": args["keywords"],
        }
        autosettings["Q"] = args["Q"]
        autosettings["pdfs"] = [original_pdfs]
        autosettings["gaussianity"] = {"pdfs": [original_pdfs]}
        autosettings["decomposition"] = {
            "pdfs": [hessian_pdfs, f"{hessian_pdfs}_pos", f"{hessian_pdfs}_neg"], "normalize_to": hessian_pdfs
        }
        autosettings["MC_Hessian_compare"] = {
            "pdfs": [hessian_pdfs, original_pdfs], "normalize_to": original_pdfs
        }

        return autosettings

 
    def get_config(self):

        def decompose(pdf_name, fit_name):

            log.info("Decomposing Hessian pdfs...")

            pdf = PDF(name=pdf_name)
            fit = fit_name
            t0pdfset = "NNPDF31_nnlo_as_0118"
            template = {
                "theory": {"from_": "fit"},
                "theoryid": {"from_": "theory"},
                "use_cuts": "fromfit",
                "experiments": {"from_": "fit"},
                "fit": fit,
                "use_t0": True,
                "t0pdfset": t0pdfset,
            }

            delta_chi2 = API.delta_chi2_hessian(pdf=pdf_name, **template)
            ind_pos = np.asarray([i for i in range(len(delta_chi2)) if delta_chi2[i] >= 0])
            ind_neg = np.asarray([i for i in range(len(delta_chi2)) if i not in ind_pos])

            ind_pos, ind_neg = ind_pos + 1, ind_neg + 1
            new_pdf_from_indexes(
                pdf=pdf, indexes=ind_pos, set_name=pdf_name + "_pos", installgrid=True, hessian=True
            )
            new_pdf_from_indexes(
                pdf=pdf, indexes=ind_neg, set_name=pdf_name + "_neg", installgrid=True, hessian=True
            )

            log.info("Finished ecomposing Hessian pdfs")

        decompose(self.args["hessian_pdfs"], self.args["original_pdfs"])

        runcard = deltachitemplates.template_path
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
