"""Compatibility shims for third-party packages that lag behind our stack (Python 3.14, SciPy >= 1.14,
setuptools >= 81). Import this module BEFORE importing triceratops.

- pytransit (a triceratops dependency) imports `scipy.integrate.trapz` and `numpy.trapz`, which
  were removed in favour of `trapezoid`.
- pytransit imports `pkg_resources.resource_filename`; `pkg_resources` was removed from setuptools 81.

Each shim only adds a missing name and never replaces an existing one.
"""

from __future__ import annotations

import importlib.resources
import sys
import types

import numpy as np
import scipy.integrate

if not hasattr(scipy.integrate, "trapz"):
    scipy.integrate.trapz = scipy.integrate.trapezoid
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

try:
    import pkg_resources  # noqa: F401
except ModuleNotFoundError:
    _pr = types.ModuleType("pkg_resources")
    _pr.resource_filename = lambda package, path: str(importlib.resources.files(package) / path)
    sys.modules["pkg_resources"] = _pr
