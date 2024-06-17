"""Test related to the seletopclusts module."""

import os
import pytest
import pytest_mock  # noqa : F401
import tempfile
import math
from pathlib import Path

from . import golden_data

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.seletopclusts import HaddockModule as SeleTopClustModule  # noqa : E501
from haddock.modules.analysis.seletopclusts import DEFAULT_CONFIG
from haddock.modules.analysis.seletopclusts.seletopclusts import (
    map_clusters_models,
    select_top_clusts_models,
    rank_clust_order,
    size_clust_order,
    sort_models,
    write_selected_models,
    )


NB_CLUSTERS = 3
MDL_RANK_INCREMENT = 5


###########################
# DEFINE FIXED INPUT DATA #
###########################
@pytest.fixture
def clustered_models() -> list[PDBFile]:
    """Set of clustered PDBfiles."""
    models: list[PDBFile] = []
    # Create 3 cluster of 7, 6 and 5 models respectively
    for clt_rank in range(1, NB_CLUSTERS + 1):
        for mdl_rank in range(clt_rank + MDL_RANK_INCREMENT - 1, 0, -1):
            pdbfile = PDBFile(
                Path(golden_data, "protprot_complex_1.pdb"),
                path=golden_data,
                )
            # Add attributes
            pdbfile.clt_rank = clt_rank
            pdbfile.clt_model_rank = mdl_rank
            models.append(pdbfile)
    return models
    

@pytest.fixture
def unclustered_models() -> list[PDBFile]:
    """Set of unclustered PDBfiles."""
    models: list[PDBFile] = []
    # Create 3 cluster of 7, 6 and 5 models respectively
    for clt_rank in range(1, NB_CLUSTERS + 1):
        for _mdl_rank in range(clt_rank + MDL_RANK_INCREMENT - 1, 0, -1):
            pdbfile = PDBFile(
                Path(golden_data, "protprot_complex_1.pdb"),
                path=golden_data,
                )
            # Do NOT add attributes
            models.append(pdbfile)
    return models


@pytest.fixture
def unranked_models() -> list[PDBFile]:
    """List of unranked models."""
    return [
        PDBFile(
            Path(golden_data, "protprot_complex_1.pdb"),
            path=golden_data,
            ),
        PDBFile(
            Path(golden_data, "protprot_complex_2.pdb"),
            path=golden_data,
            ),
        ]


@pytest.fixture
def ranked_models() -> list[PDBFile]:
    """List of ranked models."""
    models: list[PDBFile] = []
    for ind, clt_model_rank in enumerate(range(2, 0, -1), start=1):
        pdbfile = PDBFile(
            Path(golden_data, f"protprot_complex_{ind}.pdb"),
            path=golden_data,
            )
        # Add clt model rank attribute
        pdbfile.clt_model_rank = clt_model_rank
        pdbfile.clt_rank = 1
        models.append(pdbfile)
    return models


@pytest.fixture
def seletopclust():
    """Test module __init__()."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        yield SeleTopClustModule(
            order=1,
            path=Path(tmpdir),
            initial_params=DEFAULT_CONFIG,
            )


class MockPreviousIO:
    """A mocking class holding specific methods."""

    def __init__(self, models):
        self.models = models

    # In the mocked method, add the arguments that are called by the original
    #  method that is being tested
    def retrieve_models(self, individualize: bool = False):
        """Provide a set of models."""
        return self.models


###################################
# Testing of class in __init__.py #
###################################
def test_confirm_installation(seletopclust):
    """Test confirm install."""
    assert seletopclust.confirm_installation() is None


def test_init(seletopclust):
    """Test __init__ function."""
    seletopclust.__init__(
        order=42,
        path=Path("0_anything"),
        initial_params=DEFAULT_CONFIG,
        )

    # Once a module is initialized, it should have the following attributes
    assert seletopclust.path == Path("0_anything")
    assert seletopclust._origignal_config_file == DEFAULT_CONFIG
    assert type(seletopclust.params) == dict
    assert len(seletopclust.params.keys()) != 0


def test_seletopclust_run(seletopclust, mocker, clustered_models):
    """Test content of _run() function from __init__.py HaddockModule class."""
    # Mock some functions
    seletopclust.previous_io = MockPreviousIO(clustered_models)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
        )
    # run main module _run() function
    module_sucess = seletopclust.run()
    assert module_sucess is None


def test_seletopclust_neg_nb_mdls(seletopclust, mocker, clustered_models):
    """Test finish_with_error due to wrong nb models parameter value."""
    # Change parameter
    seletopclust.params["top_models"] = -1
    # Mock some functions
    seletopclust.previous_io = MockPreviousIO(clustered_models)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.finish_with_error",
        side_effect=Exception('mocked error'),
        )
    with pytest.raises(Exception) as moked_finish_with_error:
        # run main module _run() function
        seletopclust.run()
    assert moked_finish_with_error.value.__str__() == 'mocked error'


def test_seletopclust_wrong_clust_param_type(
        seletopclust,
        mocker,
        clustered_models
        ):
    """Test finish_with_error due to wrong cluster parameter type."""
    # Change parameter
    seletopclust.params["top_cluster"] = '1'
    # Mock some functions
    seletopclust.previous_io = MockPreviousIO(clustered_models)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.finish_with_error",
        side_effect=Exception('mocked error'),
        )
    with pytest.raises(Exception) as moked_finish_with_error:
        # run main module _run() function
        seletopclust.run()
    assert moked_finish_with_error.value.__str__() == 'mocked error'


def test_seletopclust_unclustered(seletopclust, mocker, unclustered_models):
    """Test finish_with_error due to unclustered data."""
    # Mock some functions
    seletopclust.previous_io = MockPreviousIO(unclustered_models)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.finish_with_error",
        side_effect=Exception('mocked error'),
        )
    with pytest.raises(Exception) as moked_finish_with_error:
        # run main module _run() function
        seletopclust.run()
    assert moked_finish_with_error.value.__str__() == 'mocked error'


#########################
# DEFINE TEST FUNCTIONS #
#########################
def test_map_clusters_models(clustered_models):
    """Test grouping of models by clusters."""
    by_clusters = map_clusters_models(clustered_models)
    # Check that we indeed have 3 clusters
    assert len(by_clusters.keys()) == 3
    # Check that clusters only contain models of the same cluster
    for clt_rank, clt_models in by_clusters.items():
        for mdl in clt_models:
            assert mdl.clt_rank == clt_rank


def test_rank_clust_order(clustered_models):
    """Test cluster ordering by rank."""
    by_clusters = map_clusters_models(clustered_models)
    cluster_rank_order = rank_clust_order(by_clusters)
    # Make sure we do not miss any cluster
    assert len(cluster_rank_order) == NB_CLUSTERS
    for clt_rank_ind, clt_rank in enumerate(cluster_rank_order, start=1):
        assert clt_rank_ind == clt_rank


def test_size_clust_order(clustered_models):
    """Test cluster ordering by size."""
    by_clusters = map_clusters_models(clustered_models)
    cluster_size_order = size_clust_order(by_clusters)
    # Make sure we do not miss any cluster
    assert len(cluster_size_order) == NB_CLUSTERS
    # Make sure the order is fine
    assert cluster_size_order[0] == 3
    assert cluster_size_order[1] == 2
    assert cluster_size_order[2] == 1


def test_select_top_ranked_clusts_models(clustered_models):
    """Test select_top_clusts_models behavior with rank ordered clusters."""
    # Test 1 model, 1 cluster
    selected_models, _notes = select_top_clusts_models(
        "score",
        clustered_models,
        1,
        1,
        )
    assert len(selected_models) == 1
    assert selected_models[0].clt_rank == 1
    assert selected_models[0].clt_model_rank == 1

    # Test max models
    selected_models, _notes = select_top_clusts_models(
        "score",
        clustered_models,
        2,
        8,
        )
    assert len(selected_models) == 5 + 6

    ind = 0
    for clt_rank in range(1, 2 + 1):
        for mdl_rank in range(1, clt_rank + MDL_RANK_INCREMENT):
            assert selected_models[ind].clt_rank == clt_rank
            assert selected_models[ind].clt_model_rank == mdl_rank
            ind += 1

    # Test max models and max clusters
    selected_models, _notes = select_top_clusts_models(
        "score",
        clustered_models,
        NB_CLUSTERS + 2,
        NB_CLUSTERS + MDL_RANK_INCREMENT + 1,
        )
    assert len(selected_models) == len(clustered_models)
    ind = 0
    for clt_rank in range(1, NB_CLUSTERS + 1):
        for mdl_rank in range(1, clt_rank + MDL_RANK_INCREMENT):
            assert selected_models[ind].clt_rank == clt_rank
            assert selected_models[ind].clt_model_rank == mdl_rank
            ind += 1

    # Test unspecified nb models, max clusters
    selected_models, _notes = select_top_clusts_models(
        "score",
        clustered_models,
        NB_CLUSTERS + 2,
        math.nan,
        )
    assert len(selected_models) == len(clustered_models)
    ind = 0
    for clt_rank in range(1, NB_CLUSTERS + 1):
        for mdl_rank in range(1, clt_rank + MDL_RANK_INCREMENT):
            assert selected_models[ind].clt_rank == clt_rank
            assert selected_models[ind].clt_model_rank == mdl_rank
            ind += 1


def test_select_top_sized_clusts_models(clustered_models):
    """Test select_top_clusts_models behavior with size ordered clusters."""
    # Test 1 model, 1 cluster
    selected_models, _notes = select_top_clusts_models(
        "size",
        clustered_models,
        1,
        1,
        )
    assert len(selected_models) == 1
    assert selected_models[0].clt_rank == 3
    assert selected_models[0].clt_model_rank == 1

    # Test max models
    selected_models, _notes = select_top_clusts_models(
        "size",
        clustered_models,
        2,
        8,
        )
    assert len(selected_models) == 7 + 6

    ind = 0
    for clt_rank in range(2, 0):
        for mdl_rank in range(1, clt_rank + MDL_RANK_INCREMENT):
            assert selected_models[ind].clt_rank == clt_rank
            assert selected_models[ind].clt_model_rank == mdl_rank
            ind += 1

    # Test max models and max clusters
    selected_models, _notes = select_top_clusts_models(
        "size",
        clustered_models,
        NB_CLUSTERS + 2,
        NB_CLUSTERS + MDL_RANK_INCREMENT + 1,
        )
    assert len(selected_models) == len(clustered_models)
    ind = 0
    for clt_rank in range(NB_CLUSTERS, 0, -1):
        for mdl_rank in range(1, clt_rank + MDL_RANK_INCREMENT):
            assert selected_models[ind].clt_rank == clt_rank
            assert selected_models[ind].clt_model_rank == mdl_rank
            ind += 1

    # Test unspecified nb models, max clusters
    selected_models, _notes = select_top_clusts_models(
        "size",
        clustered_models,
        NB_CLUSTERS + 2,
        math.nan,
        )
    assert len(selected_models) == len(clustered_models)
    ind = 0
    for clt_rank in range(NB_CLUSTERS, 0, -1):
        for mdl_rank in range(1, clt_rank + MDL_RANK_INCREMENT):
            assert selected_models[ind].clt_rank == clt_rank
            assert selected_models[ind].clt_model_rank == mdl_rank
            ind += 1


def test_sort_models_unranked(unranked_models):
    """Test sorting of unranked models."""
    sorted_models, _note = sort_models(unranked_models)
    assert sorted_models == unranked_models
    assert type(_note) == str


def test_sort_models_ranked(ranked_models):
    """Test sorting of ranked models."""
    sorted_models, _note = sort_models(ranked_models)
    assert not _note
    assert sorted_models[0].file_name == "protprot_complex_2.pdb"
    assert sorted_models[1].file_name == "protprot_complex_1.pdb"


def test_write_selected_models(ranked_models):
    """Test writing of models names mapping file."""
    with tempfile.TemporaryDirectory(dir="./") as tmpdir:
        outputfile = f"{tmpdir}test-seletopclusts.txt"
        models = write_selected_models(
            outputfile,
            ranked_models,
            './',
            )
        # Validate file creation
        assert os.path.exists(outputfile)
        assert Path(outputfile).stat().st_size != 0
        Path(outputfile).unlink(missing_ok=False)

        # Check that models were copied
        for model in models:
            outputfile = model.file_name
            assert os.path.exists(outputfile)
            assert Path(outputfile).stat().st_size != 0
            Path(outputfile).unlink(missing_ok=False)
