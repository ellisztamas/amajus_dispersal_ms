import pytest
import numpy as np

from amajusmating import mcmc

def test_check_parmeter_output():
    initial_model = {
        'missing' : 0.3, # proportion missing fathers
        'shape'  : 1,
        'scale'  : 10,
        'mixture' : 0.8
    }
    out = mcmc.check_parameters(initial_model)
    assert isinstance(out, dict)
    assert out.keys == initial_model.keys()
    assert("assortment" not in out2.keys())

    model2 = {
        'missing' : 0.3, # proportion missing fathers
        'shape'  : 1,
        'scale'  : 10,
        'mixture' : 0.8,
        'assortment' : 0.5
    }
    out2 = mcmc.check_parameters(model2)
    assert("assortment" in out2.keys())

    # Should raise ValueError because missing > 1
    model3 = {
        'missing' : 1.3, # proportion missing fathers
        'shape'  : 1,
        'scale'  : 10,
        'mixture' : 0.8,
        'assortment' : 0.5
    }
    with pytest.raises(Exception):
        mcmc.check_parameters(model3)

# Run the MCMC
# Dictionary listing starting values.
initial_model = {
    'missing' : 0.3, # proportion missing fathers
    'shape'  : 1,
    'scale'  : 30,
    'mixture' : 0.8
}