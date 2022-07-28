"""Test CNS paths creation."""
from pathlib import Path

import pytest

from haddock.core import cns_paths


@pytest.fixture(
    params=[
        cns_paths.axis,
        cns_paths.tensors,
        cns_paths.translation_vectors,
        cns_paths.water_box,
        ]
    )
def dict_containing_paths(request):
    """HADDOCK3 dict containig paths to be tested equally."""
    return request.param


@pytest.fixture(
    params=[
        # these are functions
        cns_paths.get_axis,
        cns_paths.get_tensors,
        cns_paths.get_translation_vectors,
        cns_paths.get_water_box,
        ]
    )
def dummy_generated_paths(request):
    """HADDOCK3 dict containig dummy paths to be tested equally."""
    # add the "dummy_folder" to the functions
    return request.param("dummy_folder")


def test_dict_path_is_dict(dict_containing_paths):
    """Test is dict."""
    assert isinstance(dict_containing_paths, dict)


def test_dict_path_values_are_Path(dict_containing_paths):
    """Test dict_paths are paths."""
    assert all(isinstance(p, Path) for p in dict_containing_paths.values())


def test_dict_path_values_exist(dict_containing_paths):
    """Test paths exist exist."""
    assert all(p.exists() for p in dict_containing_paths.values())


def test_dict_path_values_are_absolute_paths(dict_containing_paths):
    """Test translation vectors are absolute."""
    assert all(p.is_absolute() for p in dict_containing_paths.values())


@pytest.mark.parametrize(
    "dpaths,sw",
    [
        (cns_paths.tensors, "tensor_"),
        (cns_paths.translation_vectors, "trans_vector_"),
        ]
    )
def test_dicts_same_keys(dpaths, sw):
    """Test translation vectors startwidth keys."""
    assert all(k.startswith(sw) for k in dpaths)


@pytest.mark.parametrize(
    "dpaths,length",
    [
        (cns_paths.axis, 3),
        (cns_paths.tensors, 6),
        (cns_paths.translation_vectors, 51),
        (cns_paths.water_box, 1),
        ]
    )
def test_trans_vectors_length(dpaths, length):
    """Test translation vectors length."""
    assert len(dpaths) == length


def test_dummy_dont_exist(dummy_generated_paths):
    """Test translation vectors do not exist."""
    assert all(not p.exists() for p in dummy_generated_paths.values())


def test_dummy_relative(dummy_generated_paths):
    """Test translation vectors are relative."""
    assert all(not p.is_absolute() for p in dummy_generated_paths.values())


def test_dummy_path_parent_position(dummy_generated_paths):
    """Test translation vectors have the input path in expected place."""
    assert all(
        str(list(p.parents)[-2]) == "dummy_folder"
        for p in dummy_generated_paths.values()
        )


@pytest.mark.parametrize(
    'path',
    [
        cns_paths.link_file,
        cns_paths.scatter_lib,
        ]
    )
def test_static_file_definitions_exist(path):
    """Test if static file definitions exist."""
    assert path.exists()
