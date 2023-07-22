import pytest

from .iqtree2 import Iqtree2
from .iqtree1 import Iqtree1
from .cache import Cache
from numpy.testing import assert_allclose

def test_cache():
    
    def test_func(alignment_file, parameters):
        return {"alignment_file": alignment_file, "parameters": parameters}

    cache = Cache("test_cache.json")
    result = cache.get("file1", "param1", lambda: test_func("file1", "param1"))
    assert result == {"alignment_file": "file1", "parameters": "param1"}
    result = cache.get("file1", "param1", lambda: test_func("not file1", "not param1"))
    assert result == {"alignment_file": "file1", "parameters": "param1"}
    
@pytest.mark.parametrize("data_files", [(["example.phy"])], indirect=True)
def test_iqtree1(temp_dir, data_files):
    iqtree1_result = Iqtree1().process(temp_dir / 'data' / data_files[0], "-cmin 2 -m TEST -nt AUTO ")
    lnL = iqtree1_result.checkpoint.log_likelihood
    assert_allclose(lnL,-11224, rtol=10) 
    
@pytest.mark.parametrize("data_files", [(["example.phy"])], indirect=True)
def test_iqtree1(temp_dir, data_files):
    iqtree2_result = Iqtree2().process(temp_dir / 'data' / data_files[0], "-cmin 2 -m TEST -nt AUTO ")
    lnL = iqtree2_result.checkpoint.log_likelihood
    assert_allclose(lnL,-11224, rtol=10)     