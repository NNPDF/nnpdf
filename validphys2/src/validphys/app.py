"""
app.py

Mainloop of the validphys application.  Here we define tailoted extensions to
the reporthengine application (such as extra command line flags). Additionally
the *provider modules* that serve as source to the validphys actions are
declared here.

The entry point of the validphys application is the ``main`` funcion of this
module.
"""
import sys
import os
import logging
import contextlib


import lhapdf
from reportengine import app

from validphys.config import Config, Environment
from validphys import uploadutils
from validphys import mplstyles


providers = [
    "validphys.results",
    "validphys.commondata",
    "validphys.pdfgrids",
    "validphys.pdfplots",
    "validphys.dataplots",
    "validphys.fitdata",
    "validphys.arclength",
    "validphys.sumrules",
    "validphys.reweighting",
    "validphys.kinematics",
    "validphys.correlations",
    "validphys.chi2grids",
    "validphys.eff_exponents",
    "validphys.paramfits.dataops",
    "validphys.paramfits.plots",
    "validphys.theorycovariance.construction",
    "validphys.theorycovariance.output",
    "validphys.theorycovariance.tests",
    "validphys.replica_selector",
    "validphys.closuretest",
    "validphys.mc_gen",
    "validphys.theoryinfo",
    "validphys.pseudodata",
    "validphys.renametools",
    "validphys.covmats",
    "validphys.hyperoptplot",
    "validphys.deltachi2",
    "validphys.n3fit_data",
    "validphys.mc2hessian",
    "reportengine.report",
]

log = logging.getLogger(__name__)


class App(app.App):

    environment_class = Environment
    config_class = Config

    critical_message = """A critical error ocurred. This is likely due to one of the following reasons:

 - A bug in validphys.
 - Corruption of the provided resources (e.g. incorrect plotting files).
 - Cosmic rays hitting your CPU and altering the registers.

The traceback above should help determine the cause of the problem. If you
believe this is a bug in validphys (please discard the cosmic rays first),
please open an issue on GitHub<https://github.com/NNPDF/nnpdf/issues>,
including the contents of the following file:

%s
"""

    @property
    def default_style(self):
        return os.fspath(mplstyles.smallstyle)

    def __init__(self, name="validphys", providers=providers):
        super().__init__(name, providers)

    @property
    def argparser(self):
        parser = super().argparser

        cout = parser.add_mutually_exclusive_group()
        # We want True False or None, so that none defaults to debug or quiet.
        # That's why we use store_const
        cout.add_argument(
            "--cout",
            action="store_const",
            const=True,
            help="display C output. Default depends on log level",
        )
        cout.add_argument("--no-cout", dest="cout", action="store_const", const=False)

        net = parser.add_mutually_exclusive_group()
        net.add_argument(
            "--net",
            action="store_true",
            default=True,
            help="Enable remote loader. " "Try to download missing resources. This is the default",
        )
        net.add_argument(
            "--no-net",
            dest="net",
            action="store_false",
            help="Disable remote loader. Use only local resources.",
        )

        parser.add_argument(
            "--upload",
            action="store_true",
            help="Upload the resulting output folder to the Milan server.",
        )

        return parser

    def init(self):
        super().init()
        cout = self.args["cout"]
        if cout is None:
            if self.args["loglevel"] <= logging.DEBUG:
                cout = True
        if not cout:
            import NNPDF

            NNPDF.SetVerbosity(0)
            lhapdf.setVerbosity(0)

    @staticmethod
    def upload_context(do_upload, output):
        """If do_upload is False, do notihing. Otherwise, on enter, check the
        requiements for uploading and on exit,
        upload the output path if do_upload is True. Otherwise do nothing.
        Raise SystemExit on error."""
        if do_upload:
            return uploadutils.ReportUploader().upload_or_exit_context(output)
        return contextlib.ExitStack()

    def run(self):
        if sys.version_info < (3, 6):
            log.warning(
                "validphys 2 is discontinued on Python<3.6 and will "
                "not be longer updated. Please run\n"
                "conda install python=3.6\n\n"
                "If you have any problems, please open an issue "
                "on https://github.com/NNPDF/nnpdf/issues."
            )
        with self.upload_context(self.args["upload"], self.args["output"]):
            super().run()


def main():
    a = App()
    a.main()


if __name__ == "__main__":
    main()
