import faulthandler

import numpy as np
from hypothesis.strategies import floats, integers, sampled_from
from hypothesis import given

from validphys.core import PDF
from validphys.pdfbases import ALL_FLAVOURS

from NNPDF.ffi import lib, ffi
from NNPDF import pdf_pointer

faulthandler.enable()

@given(
    x=floats(min_value=0, max_value=1,),
    q=floats(min_value=2, max_value=1000,),
    n=integers(min_value=0, max_value=99),
    fl = sampled_from(ALL_FLAVOURS),
)
def test_fast_xfxQ_equal(x,q, n, fl):
    loaded_pdf = PDF("NNPDF31_nlo_as_0118").load()
    pointer = ffi.cast('void *', pdf_pointer(loaded_pdf))
    orig_value = loaded_pdf.xfxQ(x,q,n,fl)
    assert(lib.xfxQ(pointer, x, q, n, fl) == orig_value)
    user_data = ffi.new('struct along_x_user_data*',
            {'pdf':pointer, 'Q':q, 'member':n, 'fl':fl})
    assert lib.xfxQ_along_x(x, user_data) == orig_value
    if np.isfinite(orig_value) and x>0:
        pdfval = orig_value/x
        assert np.allclose(
            pdfval,
            lib.fxQ_along_x(x, user_data),
        )
    if -fl in ALL_FLAVOURS:
        original_xvalence = orig_value - loaded_pdf.xfxQ(x,q,n, -fl)
        xvalence = lib.xfxQ_valence_along_x(x, user_data)
        assert(xvalence == original_xvalence or
              (np.isnan(original_xvalence) and np.isnan(xvalence) ))
        if np.isfinite(xvalence) and x > 0:
            assert np.allclose(xvalence/x, lib.fxQ_valence_along_x(x,user_data))

    #Note this test has to go in the end because uset_data is
    #mutated by the C call.
    all_flavours_orig = sum(loaded_pdf.xfxQ(x,q,n,fl) for fl in ALL_FLAVOURS)
    all_flavours =  lib.xfxQ_along_x_sum_all(x, user_data)
    assert (all_flavours_orig == all_flavours or
           (np.isnan(all_flavours) and np.isnan(all_flavours_orig)))






