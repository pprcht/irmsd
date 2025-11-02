import numpy as np
from yourpkg import saxpy

def test_saxpy_basic():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    y = np.array([10.0, 20.0, 30.0, 40.0])
    out = saxpy(0.5, x, y)
    assert np.allclose(out, [10.5, 21.0, 31.5, 42.0])
    # in-place behavior
    assert out is y

