"""Unit test package for cordex."""

import importlib

import pytest


def _importorskip(modname):
    try:
        importlib.import_module(modname)
        has = True
    except ImportError:  # pragma: no cover
        has = False
    func = pytest.mark.skipif(not has, reason=f"requires {modname}")
    return has, func


has_pydruint, requires_pydruint = _importorskip("pydruint")
has_pyintorg, requires_pyintorg = _importorskip("pyintorg")
