import unittest

from NNPDF.nnpdfdb import IndexDB, map_string_string

class TestIndexDB(unittest.TestCase):
    def setUp(self):
        self.db = IndexDB('../../../nnpdfcpp/data/theory.db', "TheoryIndex")

    def test_extractmap(self):
        m = map_string_string()
        self.db.ExtractMap(1, ['PTO'], m)
        self.assertEqual(m.asdict().keys(), {'PTO'})
