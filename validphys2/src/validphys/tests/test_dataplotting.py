# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 21:27:51 2016

@author: Zahari Kassabov
"""
import unittest

from validphys.loader import Loader

from validphys.plotoptions import dataplotting

#TODO: Do something sensible with the paths here
l = Loader('../nnpdfcpp/data/', '')

config = """
x: k3
kinematics_override: dummy_transform
line_by:
  - k2

figure_by:
  - idat2bin
  - high_xq

extra_labels:
    idat2bin:  [0, 0, 0, 0, 0, 0, 0, 0, 100, 100, 100, 100, 100, 200, 200, 200, 300, 300, 300, 400, 400, 400, 500, 500, 600, 600, 700, 700, 800, 800, 900, 1000, 1000, 1100]

"""

class TestDataPlotting(unittest.TestCase):
    def test_all_plotinfo(self):
        cds = []
        for ds in l.available_datasets:
            try:
                cd = l.get_commondata(ds, 3)
                cds.append(cd)
            except:
                #We have to exclude all partially implemented sets and wrong sysnum
                continue

            dataplotting.PlotInfo.from_commondata(cd)

    def test_file(self):
        cd = l.get_commondata('HERA1CCEP', 2)
        info = dataplotting.PlotInfo.from_commondata(cd, config)
        self.assertEqual(info.extra_labels['idat2bin'][10], 100)
        tb = dataplotting.kitable(cd, info)
        #tb.columns = [*dataplotting.get_plot_kinlabels(cd), *tb.columns[3:]]
        for split_vals, data in tb.groupby(info.figure_by):
            print(', '.join('{key} = {value}'.format(key=key, value=value) for key, value in zip(info.figure_by, split_vals)))



if __name__ == '__main__':
    unittest.main()
