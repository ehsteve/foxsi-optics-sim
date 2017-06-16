"""
Created on Jul 19, 2011

@author: rtaylor
"""
from math import pi, acos, atan, tan
import numpy as np
from numpy.random import random
from foxsisim.reflectivity import Reflectivity
from scipy.integrate import quad
from scipy import stats
import astropy.units as u

mirror_reflectivity = Reflectivity()


def normalize(vector):
    """ Normalize a vector to become a unit vector."""
    vector = np.array(vector)
    if len(vector.shape) == 2:
        return vector / np.linalg.norm(vector, axis=1)[:,None]
    else:
        return vector / np.linalg.norm(vector)


def norm(vector):
    """
    Provide the norm or length of a vector. Is quantity compatible.
    """
    if len(vector.shape) == 2:
        return np.sqrt(np.sum(vector * vector, axis=1))
    else:
        return np.sqrt(np.sum(vector * vector))


def angle_between(a, b):
    """
    Returns the angle between two vectors. Vector can have or not have dimensions.
    """
    a_u = normalize(a)
    b_u = normalize(b)
    if len(a_u.shape) == 2:
        return np.rad2deg(np.arccos(np.clip(np.sum(a_u * b_u, axis=1), -1.0, 1.0))) * u.deg
    else:
        return np.rad2deg(np.arccos(np.clip(np.sum(a_u * b_u), -1.0, 1.0))) * u.deg


def grazing_angle(x, normal):
    """
    Calculate the grazing angle between a vector and the surface.
    The complement to the angle of incidence.
    """
    return 90 * u.deg - angle_between(x, normal)


def reflect(x, normal, energy=None):
    """
    Takes a vector and reflects it in relation to a normal
    """
    # if the angle of incidence is greater than 90 deg, return None
    incident_angle = angle_between(x, normal)
    #graze_angle = 90 * u.deg - incident_angle
    #if incident_angle > 90 * u.deg:
    #    return None
    # TODO add removal of rays if hitting back of shell
    if energy is not None:
        if random(1)[0] > mirror_reflectivity.value(energy, np.rad2deg(graze_angle)):
            return None
    if len(x.shape) == 2:
        result = normalize(x - 2. * np.sum(x * normal, axis=1)[:, None] * normal)
    return result


@u.quantity_input(radius=u.cm, focalLength=u.cm)
def calcShellAngle(radius, focal_length):
    """
    Calculates shell angle(s) given one or more shell radius
    and a focal length. Takes one radius, a list, or a numpy array
    of radii. Returns an angle or a numpy array of angles.
    """
    # took this eqn from the excel file...
    return (0.25 * np.arctan(radius / focal_length)).to(u.deg)


@u.quantity_input(radius=u.deg, focalLength=u.cm)
def calcShellRadius(angle, focal_length):
    """
    Calculates shell radii given one or more shell angles
    and a focal length. Takes one angle, a list, or a numpy array
    of angles. Returns a radius or a numpy array of radii.
    """
    return np.tan(4 * angle) * focal_length


def genCustomRands(f, n):
    """
    Generate n random numbers given a distribution (x, y)
    """
    dist = getDistribution(f)
    return dist.rvs(size=n)


def getDistribution(f):
    """
    Create a distribution to draw random numbers from
    """
    class rv(stats.rv_continuous):
        def _pdf(self, x):
            return f(x)

        def _cdf(self, x):
            return quad(f, 0, x)

    return rv(name='customdist')


def to_cm(list_or_array):
    return u.Quantity(list_or_array).to('cm').value
