#!/usr/bin/env python3
"""Check source timing, wall accounting and external fields in the coupled run."""

import sys
from pathlib import Path

import numpy as np
import scipy.constants as sc
import yt

yt.funcs.mylog.setLevel(50)
max_level = int(sys.argv[1])
charge = 100 * sc.e


def particle_count(data, species):
    if (species, "particle_weight") not in data.ds.field_list:
        return 0
    return data[species, "particle_weight"].size


for step in range(5):
    ds = yt.load(str(Path("diags") / f"diag{step:06d}"))
    assert ds.index.max_level == max_level
    data = ds.all_data()
    # At t=0 these counts catch injection after the initial field solve.
    # A deliberately mismatched my_constants.dt catches stale injection dt.
    assert particle_count(data, "electrons") == 2 + int(2.3 * step)
    assert particle_count(data, "ions") == (2 if step == 0 else 0)
    assert particle_count(data, "neutrals") == (2 if step == 0 else 4)
    for component in ("Ex", "Ey", "Ez", "phi"):
        values = data["boxlib", component].to_ndarray()
        assert np.all(np.isfinite(values))
    if step == 0:
        # Separated initial electron/ion sources must already produce a field.
        assert np.max(np.abs(data["boxlib", "Ez"].to_ndarray())) > 1.0e-6
    for grid in ds.index.grids:
        # Check each level away from PEC faces (which impose their own field BC).
        interior = np.ones(grid.ActiveDimensions, dtype=bool)
        for axis in ("x", "y", "z"):
            interior &= np.abs(grid["index", axis].to_ndarray()) < 0.0035
        if np.any(interior):
            np.testing.assert_allclose(
                grid["boxlib", "Bz"].to_ndarray()[interior],
                0.001,
                rtol=1.0e-6,
                atol=1.0e-12,
            )

wall = np.loadtxt("wall_charge.dat", ndmin=2)
assert wall.shape == (4, 7 + max_level)
np.testing.assert_array_equal(wall[:, 0], np.arange(1, 5))
np.testing.assert_allclose(wall[:, 1], np.arange(1, 5) * 1.0e-10, rtol=1.0e-12)
# Two wall impacts, each contributing +100e, with no repeated accumulation.
np.testing.assert_allclose(wall[:, 2:4], 2 * charge, rtol=1.0e-6, atol=1.0e-25)
expected_per_level = 2 * charge if max_level == 0 else charge
np.testing.assert_allclose(
    wall[:, 4 : 5 + max_level], expected_per_level, rtol=1.0e-6, atol=1.0e-25
)
# Impacts are deliberately far from refinement interfaces and physical faces.
np.testing.assert_allclose(wall[:, -2], 2 * charge, rtol=1.0e-6, atol=1.0e-25)
np.testing.assert_allclose(wall[:, -1], 0, atol=2 * charge * 1.0e-6)
print("Hall integration: initial sources, injection dt, wall charge and B passed.")
