#!/usr/bin/env python3

import sys

import numpy as np
import scipy.constants as scc
import yt


def make_axis_frame(axis):
    z_axis = np.asarray(axis, dtype=float)
    z_axis /= np.linalg.norm(z_axis)
    denominator = 1.0 + z_axis[2]
    if denominator <= 16.0 * np.finfo(float).eps:
        return np.array(
            [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]
        ).T

    inverse_denominator = 1.0 / denominator
    x_axis = np.array(
        [
            1.0 - z_axis[0] ** 2 * inverse_denominator,
            -z_axis[0] * z_axis[1] * inverse_denominator,
            -z_axis[0],
        ]
    )
    y_axis = np.array(
        [
            -z_axis[0] * z_axis[1] * inverse_denominator,
            1.0 - z_axis[1] ** 2 * inverse_denominator,
            -z_axis[1],
        ]
    )
    return np.column_stack((x_axis, y_axis, z_axis))


dataset = yt.load(sys.argv[1])
data = dataset.all_data()

x = data["neutrals", "particle_position_x"].to_ndarray()
y = data["neutrals", "particle_position_y"].to_ndarray()
px = data["neutrals", "particle_momentum_x"].to_ndarray()
py = data["neutrals", "particle_momentum_y"].to_ndarray()
pz = data["neutrals", "particle_momentum_z"].to_ndarray()

x_relative = x - 0.35
y_relative = y - 0.45
radius = np.sqrt(x_relative**2 + y_relative**2)
theta = np.arctan2(y_relative, x_relative)

proper_velocity = np.vstack((px, py, pz)) / scc.m_p
gamma = np.sqrt(1.0 + np.sum(proper_velocity**2, axis=0) / scc.c**2)
velocity = proper_velocity / gamma

cos_theta = np.cos(theta)
sin_theta = np.sin(theta)
velocity_at_theta0 = np.vstack(
    (
        velocity[0] * cos_theta + velocity[1] * sin_theta,
        -velocity[0] * sin_theta + velocity[1] * cos_theta,
        velocity[2],
    )
)
local_velocity = make_axis_frame([3.0, 4.0, 12.0]).T @ velocity_at_theta0

assert x.size == 20000
assert np.all((radius >= 0.1) & (radius <= 0.2))
assert np.all(local_velocity[2] >= 0.0)

expected_means = np.array([0.0, 0.0, 2.0e5])
expected_sigmas = np.array([1.0e4, 2.0e4, 3.0e4])
measured_means = np.mean(local_velocity, axis=1)
measured_sigmas = np.std(local_velocity, axis=1)

mean_tolerances = np.array([5.0e2, 1.0e3, 1.0e3])
assert np.all(np.abs(measured_means - expected_means) < mean_tolerances)
assert np.all(np.abs(measured_sigmas / expected_sigmas - 1.0) < 0.03)
