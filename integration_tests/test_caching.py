"""Integration test for the caching functionality."""

import pytest
import tempfile
from pathlib import Path
import pickle

from haddock.libs.libontology import ModuleIO, Cache, CacheKey


def test_cache_basic_operations():
    """Test basic cache operations: add, get, has, get_all."""
    cache = Cache()

    # Test adding items
    cache.add(CacheKey.RMSD, ("model1", "model2"), 1.5)
    cache.add(CacheKey.RMSD, ("model1", "model3"), 2.3)
    cache.add(CacheKey.RMSD, ("model2", "model3"), 0.8)

    # Test getting items
    assert cache.get(CacheKey.RMSD, ("model1", "model2")) == 1.5
    assert cache.get(CacheKey.RMSD, ("model1", "model3")) == 2.3
    assert cache.get(CacheKey.RMSD, ("model2", "model3")) == 0.8

    # Test has method
    assert cache.has(CacheKey.RMSD, ("model1", "model2"))
    assert not cache.has(CacheKey.RMSD, ("model3", "model4"))

    # Test get_all
    all_rmsd = cache.get_all(CacheKey.RMSD)
    assert len(all_rmsd) == 3


def test_cache_clear():
    """Test clearing cache entries."""
    cache = Cache()

    # Add some items
    cache.add(CacheKey.RMSD, ("model1", "model2"), 1.5)
    cache.add(CacheKey.RMSD, ("model3", "model4"), 2.0)

    # Verify items exist
    assert cache.has(CacheKey.RMSD, ("model1", "model2"))
    assert cache.has(CacheKey.RMSD, ("model3", "model4"))
    assert len(cache.get_all(CacheKey.RMSD)) == 2

    # Clear RMSD cache
    cache.clear(CacheKey.RMSD)

    # Verify cache is empty
    assert not cache.has(CacheKey.RMSD, ("model1", "model2"))
    assert not cache.has(CacheKey.RMSD, ("model3", "model4"))
    assert len(cache.get_all(CacheKey.RMSD)) == 0


def test_module_io_cache_serialization():
    """Test that ModuleIO can save and load cache data."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create ModuleIO instance and add some cache data
        module_io = ModuleIO()
        module_io.cache.add(CacheKey.RMSD, ("model1", "model2"), 1.5)
        module_io.cache.add(CacheKey.RMSD, ("model3", "model4"), 2.0)

        # Save to file
        save_path = Path(tmpdir) / "test_io.json"
        module_io.save(save_path.parent, save_path.name)

        # Verify file was created
        assert save_path.exists()
        assert save_path.stat().st_size > 0

        # Create new ModuleIO instance and load
        new_module_io = ModuleIO()
        new_module_io.load(save_path)

        # Verify cache was loaded correctly
        assert new_module_io.cache.has(CacheKey.RMSD, ("model1", "model2"))
        assert new_module_io.cache.get(CacheKey.RMSD, ("model1", "model2")) == 1.5
        assert new_module_io.cache.has(CacheKey.RMSD, ("model3", "model4"))
        assert new_module_io.cache.get(CacheKey.RMSD, ("model3", "model4")) == 2.0

        # Verify all cache entries
        all_rmsd = new_module_io.cache.get_all(CacheKey.RMSD)
        assert len(all_rmsd) == 2


def test_cache_pickle_serialization():
    """Test that Cache can be pickled and unpickled directly."""
    cache = Cache()

    # Add test data
    cache.add(CacheKey.RMSD, ("model1", "model2"), 1.5)
    cache.add(CacheKey.RMSD, ("model3", "model4"), 2.0)

    # Pickle the cache
    pickled_data = pickle.dumps(cache)

    # Unpickle and verify
    new_cache = pickle.loads(pickled_data)

    assert new_cache.has(CacheKey.RMSD, ("model1", "model2"))
    assert new_cache.get(CacheKey.RMSD, ("model1", "model2")) == 1.5
    assert new_cache.has(CacheKey.RMSD, ("model3", "model4"))
    assert new_cache.get(CacheKey.RMSD, ("model3", "model4")) == 2.0

    # Verify all entries
    all_rmsd = new_cache.get_all(CacheKey.RMSD)
    assert len(all_rmsd) == 2


def test_cache_multiple_keys():
    """Test cache with multiple CacheKey types."""
    cache = Cache()

    # Add RMSD entries
    cache.add(CacheKey.RMSD, ("model1", "model2"), 1.5)
    cache.add(CacheKey.RMSD, ("model3", "model4"), 2.0)

    # Verify RMSD entries don't interfere with each other
    rmsd_entries = cache.get_all(CacheKey.RMSD)
    assert len(rmsd_entries) == 2

    # Test that we can clear specific cache keys independently
    cache.clear(CacheKey.RMSD)
    assert len(cache.get_all(CacheKey.RMSD)) == 0


def test_cache_nonexistent_key():
    """Test cache behavior with non-existent keys."""
    cache = Cache()

    # Test getting non-existent key
    result = cache.get(CacheKey.RMSD, ("nonexistent1", "nonexistent2"))
    assert result is None

    # Test has with non-existent key
    assert not cache.has(CacheKey.RMSD, ("nonexistent1", "nonexistent2"))

    # Test get_all with empty cache
    all_rmsd = cache.get_all(CacheKey.RMSD)
    assert len(all_rmsd) == 0


def test_cache_update_method():
    """Test the update_cache method in ModuleIO."""
    # Create two ModuleIO instances
    module_io1 = ModuleIO()
    module_io2 = ModuleIO()

    # Add cache data to first instance
    module_io1.cache.add(CacheKey.RMSD, ("model1", "model2"), 1.5)

    # Update second instance with first instance's cache
    module_io2.update_cache(module_io1.cache)

    # Verify the update worked
    assert module_io2.cache.has(CacheKey.RMSD, ("model1", "model2"))
    assert module_io2.cache.get(CacheKey.RMSD, ("model1", "model2")) == 1.5


def test_cache_with_module_io_lifecycle():
    """Test cache behavior through complete ModuleIO lifecycle."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create and populate ModuleIO
        original_io = ModuleIO()
        original_io.cache.add(CacheKey.RMSD, ("model1", "model2"), 1.5)
        original_io.cache.add(CacheKey.RMSD, ("model3", "model4"), 2.0)

        # Add some input/output files (to test they don't interfere)
        from haddock.libs.libontology import PDBFile

        test_pdb = Path(tmpdir) / "test.pdb"
        test_pdb.write_text("ATOM      1  N   ALA A   1       1.000   1.000   1.000")
        original_io.add(PDBFile(test_pdb, path=tmpdir))

        # Save and load
        save_path = Path(tmpdir) / "test_io.pkl"
        original_io.save(save_path.parent, save_path.name)

        loaded_io = ModuleIO()
        loaded_io.load(save_path)

        # Verify cache survived the round trip
        assert loaded_io.cache.has(CacheKey.RMSD, ("model1", "model2"))
        assert loaded_io.cache.get(CacheKey.RMSD, ("model1", "model2")) == 1.5
        assert loaded_io.cache.has(CacheKey.RMSD, ("model3", "model4"))
        assert loaded_io.cache.get(CacheKey.RMSD, ("model3", "model4")) == 2.0

        # Verify input/output also survived
        assert len(loaded_io.input) == 1
        assert len(loaded_io.output) == 0
