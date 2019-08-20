import pytest
from preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
from preprocessor.meshHandle.dualCoarseMesh import DualCoarseMesh as dual
import preprocessor.geoUtil.geoTools as gtool

@pytest.fixture
def get_elements():
    M = msh('mesh/test_mesh10.h5m', dim = 3)
    volumes = M.volumes
    return volumes

def test_finescale_elements(get_elements):
    volumes = get_elements()
    assert len(volumes) == 1000
