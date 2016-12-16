# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 10:19:35 2016

@author: Zahari Kassabov
"""
import pathlib
import logging

from reportengine import app

from validphys.config import Config, Environment
#from validphys import providers


providers = [
             'validphys.results',
             'validphys.plots',
             'validphys.reweighting',
             'validphys.fitdata',
             'validphys.pdfgrids',
             'reportengine.report'
            ]

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

        parser.add_argument('-p','--datapath', help="path where the NNPDF "
                        "data is located",
                        default='../nnpdfcpp/data')

        parser.add_argument('--resultspath', help="path where the fit results "
                          "are located. Calculated from 'datapath' by default",
                         )

        cout = parser.add_mutually_exclusive_group()
        #We want True False or None, so that none defaults to debug or quiet.
        #That's why we use store_const
        cout.add_argument('--cout', action='store_const', const=True,
                          help = "display C output. Default depends on log level")
        cout.add_argument('--no-cout', dest='cout',
                              action='store_const', const=False)

        return parser



    def get_commandline_arguments(self):
        args = super().get_commandline_arguments()
        if not args['resultspath']:
            args['resultspath'] = pathlib.Path(args['datapath']).parent / 'nnpdfbuild' / 'results'
        return args

    def init(self):
        super().init()
        cout = self.args['cout']
        if cout is None:
            if self.args['loglevel'] <= logging.DEBUG:
                cout = True
        if not cout:
            from NNPDF import common
            from NNPDF.lhapdfset import setVerbosity
            common.SetVerbosity(0)
            #No idea why this doesn't work
            #import lhapdf
            setVerbosity(0)

def main():
    a = App()
    a.init()
    a.run()

if __name__ == '__main__':
    main()