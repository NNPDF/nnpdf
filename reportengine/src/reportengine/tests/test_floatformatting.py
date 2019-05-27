import decimal

from hypothesis import given
from hypothesis.strategies import (floats, lists, integers, tuples, decimals,
                                   one_of, none)
import pandas as pd
import numpy as np

from reportengine.floatformatting import (format_number, significant_digits,
                                          format_error_value_columns,
                                          ValueErrorTuple,
                                          write_in_adequate_representation)

@given(floats(allow_nan=False))
def test_format_rountrip(x):
    assert (significant_digits(x, 4) == decimal.Decimal(format_number(x, 4)))

def test_nan_roundtrip():
    assert decimal.Decimal(format_number(float('nan'))).is_nan()

samelists = integers(min_value=1, max_value=10).flatmap(lambda n:
    tuples(lists(floats(), min_size=n, max_size=n),
           lists(floats(), min_size=n, max_size=n))
)
@given(samelists)
def test_cols(valerr):
    df = pd.DataFrame(np.array(valerr).T, columns=['Val','Err'])
    formatted = format_error_value_columns(df, 'Val', 'Err')
    format_error_value_columns(df, 'Val', 'Err', inplace=True)
    assert (formatted == df).all().all()

@given(floats(), floats())
def test_valueerrortuple(value, error):
    assert '±' in str(ValueErrorTuple(value,error))

int_or_none = one_of(integers(), none())

@given(decimals(allow_nan=False, allow_infinity=False),int_or_none, int_or_none)
def test_correct_writing(d, minexp, maxexp):
    assert decimal.Decimal(d) == decimal.Decimal(write_in_adequate_representation(d, minexp=minexp, maxexp=maxexp))

def test_can_format_numpy():
    x = np.int64(6)
    assert format_number(x) == '6'
