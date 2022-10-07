import pytest
import numpy as np

exec(open('03_analysis/01_data_formatting/setup_FAPS_GPS.py').read())

def test_there_are_60_mothers():
    assert len(am_data.paternity) == 60

def test_missing_parents_update():
    assert am_data.paternity['J1246'].missing_parents == 0
    am_data.update_missing_dads(0.1)
    assert am_data.paternity['J1246'].missing_parents == 0.1
    with pytest.raises(ValueError):
        am_data.update_missing_dads(x=1.1)

def test_update_dispersal():
    assert am_data.covariates == {}
    assert am_data.paternity['J1246'].covariate == 0

    am_data.update_dispersal_probs(scale=10.0, shape=2.0, mixture = 0.8)
    assert am_data.covariates['dispersal'].shape == (60,2094)
    assert np.isfinite(am_data.covariates['dispersal']).all()

    am_data.update_dispersal_probs(scale=10.0, shape=2.0, mixture = 0.8, max_distance= 500)
    assert np.isinf(am_data.covariates['dispersal']).any()
    assert not np.isinf(am_data.covariates['dispersal']).all()

def test_update_assortment():
    am_data.update_assortment_probs(0.1)
    assert am_data.covariates['assortment'].shape == (60,2094)
    assert np.isfinite(am_data.covariates['assortment'])

def test_sibship_clustering():
    with pytest.raises(KeyError):
        am_data.params['loglik']
    
    am_data.sibship_clustering(use_covariates=True)
    assert isinstance(am_data.sibships, dict)
    assert am_data.sibships.keys() == am_data.mothers
    assert 'loglik' in am_data.params
    assert am_data.params['loglik'] < 0

    