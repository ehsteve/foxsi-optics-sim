import pytest
import foxsisim.mymath as m
from numpy.random import uniform
from numpy.testing import assert_allclose
import astropy.units as u
import numpy as np


def test_normalizing_vector_vectorized():
    dim = 100
    random_vectors = uniform(0, 1, (dim, 3)) * u.cm
    norm_random_vectors = m.unit_vector(random_vectors)
    assert_allclose(np.linalg.norm(norm_random_vectors, axis=1), np.ones(dim))


def test_normalizing_vector_nonvectorized():
    dim = 100
    random_vectors = uniform(0, 1, (dim, 3)) * u.cm
    norm_random_vectors = m.unit_vector(random_vectors)
    for this_vector in norm_random_vectors:
        assert_allclose(np.linalg.norm(this_vector), 1)
    assert_allclose(np.linalg.norm(norm_random_vectors, axis=1), np.ones(dim))


def test_angle_between_nonvectorized():
    assert_allclose(m.angleBetweenVectors((1, 0, 0)*u.cm, (0, 1, 0)*u.cm), 1.5707963267948966*u.rad)
    assert_allclose(m.angleBetweenVectors((1, 0, 0)*u.cm, (1, 0, 0)*u.cm), 0.0*u.rad)
    assert_allclose(m.angleBetweenVectors((1, 0, 0)*u.cm, (-1, 0, 0)*u.cm), np.pi*u.rad)

