"""
test_vplistscript.py

Module for testing vp-list. The output of which is dynamic and so we just check
that the script runs and gives some output
"""
from contextlib import redirect_stdout
import io

from validphys.scripts.vp_list import main


def test_listfits():
    """Checks listing fits returns output"""
    f = io.StringIO()
    cmd = ["fits"]
    with redirect_stdout(f):
        main(cmd)
    assert f.getvalue()


def test_listpdfs():
    """Checks listing pdfs returns output"""
    f = io.StringIO()
    cmd = ["pdfs"]
    with redirect_stdout(f):
        main(cmd)
    assert f.getvalue()


def test_listtheories():
    """Checks listing theories returns output"""
    f = io.StringIO()
    cmd = ["theories"]
    with redirect_stdout(f):
        main(cmd)
    assert f.getvalue()


def test_listdatasets():
    """Checks listing datasets returns output"""
    f = io.StringIO()
    cmd = ["datasets"]
    with redirect_stdout(f):
        main(cmd)
    assert f.getvalue()


def test_local():
    """Check local flag"""
    f = io.StringIO()
    cmd = ["datasets", "-l"]
    with redirect_stdout(f):
        main(cmd)
    assert f.getvalue()


def test_remote():
    """Test remote flag on both datasets (which should return empty string) and pdfs
    which returns output

    """
    f = io.StringIO()
    cmd_data = ["datasets", "-r"]
    cmd_pdfs = ["pdfs", "-r"]
    with redirect_stdout(f):
        main(cmd_data)
        assert not f.getvalue()
        main(cmd_pdfs)
    assert f.getvalue()
