import pytest
import importlib.util
from visidata import VisiData, Path, vd, options
from pathlib import Path as PathLib


@pytest.fixture
def vd():
    """Provide a fresh VisiData instance for each test"""
    vd_inst = VisiData()
    vd_inst.options = options
    return vd_inst


@pytest.fixture
def test_data_dir():
    """Return path to test data directory"""
    return PathLib(__file__).parent / "data"


@pytest.fixture
def load_plugin(vd):
    """Fixture to load a plugin into VisiData instance"""

    def _load_plugin(plugin_name):
        plugin_path = (
            PathLib(__file__).parent.parent / "vd_bioinformatics" / f"{plugin_name}.py"
        )

        # Load the plugin module directly using importlib
        spec = importlib.util.spec_from_file_location(plugin_name, plugin_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Initialize the plugin if it has an init function
        init_func = getattr(module, f"init_{plugin_name}", None)
        if init_func:
            init_func(vd)

    return _load_plugin
