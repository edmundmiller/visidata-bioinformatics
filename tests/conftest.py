import pytest
from visidata import VisiData, Path
from pathlib import Path as PathLib


@pytest.fixture
def vd():
    """Provide a fresh VisiData instance for each test"""
    return VisiData()


@pytest.fixture
def test_data_dir():
    """Return path to test data directory"""
    return PathLib(__file__).parent / "data"


@pytest.fixture
def load_plugin(vd):
    """Fixture to load a plugin into VisiData instance"""

    def _load_plugin(plugin_name):
        plugin_path = PathLib(__file__).parent.parent / "plugins" / f"{plugin_name}.py"
        vd.loadPlugin(plugin_path)

    return _load_plugin
