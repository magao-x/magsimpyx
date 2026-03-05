import pytest
import numpy as np
from magsimpyx import make_magaox_bump_mask, make_magaox_large_lyot_stop, make_gmt_lyot_aperture
from hcipy import make_pupil_grid

@pytest.fixture
def pupil_grid_1():
    return make_pupil_grid(256, 1.0)

@pytest.fixture
def pupil_grid_6p5m():
    return make_pupil_grid(256, 6.5)

@pytest.fixture
def pupil_grid_9mm():
    return make_pupil_grid(256, 9e-3)    

@pytest.fixture
def pupil_grid_25p448m():
    return make_pupil_grid(256, 25.448)

def test_make_magaox_bump_mask(pupil_grid_1, pupil_grid_6p5m, pupil_grid_9mm):
    '''Test that the bump mask is correctly scaled for different pupil grid sizes.
    '''

    bump_mask_1 = make_magaox_bump_mask(normalized=True)(pupil_grid_1)
    bump_mask_6p5m = make_magaox_bump_mask()(pupil_grid_6p5m)
    bump_mask_9mm = make_magaox_bump_mask(pupil_diameter=9e-3)(pupil_grid_9mm)

    assert np.all(bump_mask_1 == bump_mask_6p5m)
    assert np.all(bump_mask_1 == bump_mask_9mm)

def test_make_magaox_large_lyot_stop(pupil_grid_1, pupil_grid_6p5m, pupil_grid_9mm):
    '''Test that the large Lyot stop is correctly scaled for different pupil grid sizes.
    '''

    lyot_stop_1 = make_magaox_large_lyot_stop(normalized=True)(pupil_grid_1)
    lyot_stop_6p5m = make_magaox_large_lyot_stop()(pupil_grid_6p5m)
    lyot_stop_9mm  = make_magaox_large_lyot_stop(pupil_diameter=9e-3)(pupil_grid_9mm)

    assert np.all(lyot_stop_1 == lyot_stop_6p5m)
    assert np.all(lyot_stop_1 == lyot_stop_9mm)

def test_make_magaox_gmt_lyot_aperture(pupil_grid_1, pupil_grid_25p448m):
    '''Test that the GMT Lyot aperture is correctly scaled for different pupil grid sizes.
    '''

    lyot_aperture_1 = make_gmt_lyot_aperture(normalized=True)(pupil_grid_1)
    lyot_aperture_25p448m = make_gmt_lyot_aperture()(pupil_grid_25p448m)

    assert np.all(lyot_aperture_1 == lyot_aperture_25p448m)
