import pytest
import foxsisim.mymath as m
from numpy.random import uniform
from numpy.testing import assert_allclose
import astropy.units as u
import numpy as np


def test_normalizing_vector_vectorized():
    dim = 100
    random_vectors = uniform(0, 1, (dim, 3)) * u.cm
    norm_random_vectors = m.normalize(random_vectors)
    assert_allclose(np.linalg.norm(norm_random_vectors, axis=1), np.ones(dim))


def test_normalizing_vector_nonvectorized():
    dim = 100
    random_vectors = uniform(0, 1, (dim, 3)) * u.cm
    norm_random_vectors = m.normalize(random_vectors)
    for this_vector in norm_random_vectors:
        assert_allclose(np.linalg.norm(this_vector), 1)
    assert_allclose(np.linalg.norm(norm_random_vectors, axis=1), np.ones(dim))


@pytest.mark.parametrize('a, b, result', [
    ((1, 0, 0), (0, 1, 0), 0.5 * np.pi * u.rad),
    ((1, 0, 0), (1, 0, 0), 0.0 * u.rad),
    ((1, 0, 0), (-1, 0, 0), np.pi * u.rad),
])
def test_angle_between_nonvectorized(a, b, result):
    assert_allclose(m.angle_between(a,b), result)


def test_angle_between_vectorized():
    a = [(1, 0, 0), (1, 0, 0), (1, 0, 0)] * u.cm
    b = [(0, 1, 0), (1, 0, 0), (-1, 0, 0)] * u.cm
    assert_allclose(m.angle_between(a, b), [0.5 * np.pi, 0.0, np.pi] * u.rad)


@pytest.mark.parametrize('radius, focal_length, angle', [
    (1*u.cm, 1*u.cm, 0.25 * 45 * u.deg),
    (0 * u.cm, 0.00001 * u.cm, 0.0 * u.deg),
    (1*u.m, 10*u.km, 0.25 * 1/10000. * u.rad),
])
def test_calcShellAngle(radius, focal_length, angle):
    assert_allclose(m.calcShellAngle(radius, focal_length), angle)


@pytest.mark.parametrize('radius,focal_length,angle', [
    (1*u.cm, 1*u.cm, 0.25 * 45 * u.deg),
    (0 * u.cm, 0.00001 * u.cm, 0.0 * u.deg),
    (1*u.m, 10*u.km, 0.25 * 1/10000. * u.rad),
])
def test_calcShellRadius(radius, focal_length, angle):
    assert_allclose(m.calcShellRadius(angle, focal_length), radius)


@pytest.mark.parametrize('vector, result', [
    ((1, 1, 1) * u.cm, np.sqrt(3) * u.cm),   # with units
    (np.array((1, 1, 1,)), np.sqrt(3)),        # with no units
    ((3, 4, 0) * u.cm, 5 * u.cm)
])
def test_norm(vector, result):
    assert_allclose(m.norm(vector), result)


def test_norm_vectorized():
    dim = 100
    random_vectors = uniform(0, 1, (dim, 3)) * u.cm
    assert(len(m.norm(random_vectors)), dim)


@pytest.mark.parametrize('input, result', [
    ((1, 1, 1) * u.cm, np.array((1,1,1))),   # with vector
    (4 * u.cm, 4)
])
def test_to_cm(input, result):
    assert_allclose(m.to_cm(input), result)