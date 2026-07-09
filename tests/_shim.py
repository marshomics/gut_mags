"""Test shim.

The pipeline stores intermediate tables as parquet. Where pyarrow is not
installed (e.g. a minimal CI container) the tests transparently back parquet with
pickle so the analysis code runs unchanged. On the cluster, pyarrow is present
and this shim does nothing.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "python"))

try:
    import pyarrow  # noqa: F401
    HAVE_PYARROW = True
except Exception:
    HAVE_PYARROW = False
    _read_pickle = pd.read_pickle
    pd.DataFrame.to_parquet = lambda self, path, index=False, **k: self.to_pickle(path)
    pd.read_parquet = lambda path, **k: _read_pickle(path)


def run(module_name, argv):
    """Invoke a pipeline script's main() with the given CLI args."""
    import importlib
    sys.argv = [module_name] + argv
    m = importlib.import_module(module_name)
    importlib.reload(m)
    m.main()
