# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 10:19:35 2016

@author: Zahari Kassabov
"""
import sys
import logging
import contextlib


from reportengine import app

from validphys.config import Config, Environment
from validphys import uploadutils
#from validphys import providers


providers = [
             'validphys.results',
             'validphys.plots',
             'validphys.sumrules',
             'validphys.reweighting',
             'validphys.fitdata',
             'validphys.pdfgrids',
             'validphys.kinematics',
             'validphys.correlations',
             'validphys.paramfits',
             'reportengine.report'
            ]

log = logging.getLogger(__name__)

class App(app.App):

    environment_class = Environment
    config_class = Config

    critical_message = (
"""A critical error ocurred. This is likely due to one of the following reasons:

 - A bug in validphys.
 - Corruption of the provided resources (e.g. incorrect plotting files).
 - Cosmic rays hitting your CPU and altering the registers.

The traceback above should help determine the cause of the problem. If you
believe this is a bug in validphys (please discard the cosmic rays first),
please send an email to Zahari<kassabov@to.infn.it>, including the following
file in attachment:

%s
"""
    )

    @property
    def default_style(self):
        return str(self.this_folder() / 'small.mplstyle')

    def __init__(self):
        super().__init__('validphys', providers)

    @property
    def argparser(self):
        parser = super().argparser

        parser.add_argument('-p','--datapath', help="DEPRECATED USE nnprofile "
            "INSTEAD. path where the NNPDF "
            "data is located",
            default=None)

        parser.add_argument('--resultspath', help="DEPRECATED USE nnprofile "
                         "INSTEAD. path where the fit results "
                          "are located. Calculated from 'datapath' by default",
                         )

        cout = parser.add_mutually_exclusive_group()
        #We want True False or None, so that none defaults to debug or quiet.
        #That's why we use store_const
        cout.add_argument('--cout', action='store_const', const=True,
                          help = "display C output. Default depends on log level")
        cout.add_argument('--no-cout', dest='cout',
                              action='store_const', const=False)

        net = parser.add_mutually_exclusive_group()
        net.add_argument('--net', action='store_true', default=True,
                         help="Enable remote loader. "
                         "Try to download missing resources. This is the default")
        net.add_argument('--no-net', dest='net', action='store_false',
                         help="Disable remote loader. Use only local resources.")

        parser.add_argument('--upload', action='store_true',
                            help="Upload the resulting output folder to the Milan server.")

        return parser


    def init(self):
        super().init()
        cout = self.args['cout']
        if cout is None:
            if self.args['loglevel'] <= logging.DEBUG:
                cout = True
        if not cout:
            import NNPDF
            NNPDF.SetVerbosity(0)
            #No idea why this doesn't work
            #import lhapdf
            NNPDF.setVerbosity(0)

    @staticmethod
    def upload_context(do_upload, output):
        """If do_upload is False, do notihing. Otherwise, on enter, check the
        requiements for uploading and on exit,
        upload the output path if do_upload is True. Otherwise do nothing.
        Raise SystemExit on error."""
        if do_upload:
            return uploadutils.upload_or_exit_context(output)
        return contextlib.ExitStack()

    def run(self):
        if sys.version_info < (3, 6):
            log.warn("validphys 2 is discontinued on Python<3.6 and will not be longer updated. Please run\n"
                     "conda install python=3.6\n\n If you have any problems, please ask Zahari.")
        with self.upload_context(self.args['upload'], self.args['output']):
            super().run()

def main():
    a = App()
    a.main()

if __name__ == '__main__':
    main()
