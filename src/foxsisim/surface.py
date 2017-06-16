"""
Created on Jul 8, 2011

@author: rtaylor
"""
import numpy as np
import astropy.units as u
from foxsisim.mymath import normalize

class Surface:
    """
    A surface in 3-space defined by a set of parametric equations. Inherited by
    other classes like plane and segment. Provides methods for calculating points
    and normals at parameters u and v.
    """

    def __init__(self):
        """
        Constructor
        """
        pass

    def inRange(self,u,v):
        """
        Are the u and v parameters in the desired range?
        """
        return True

    def existsInOctant(self,octant):
        """
        Returns whether the surface exists in the supplied octant ***OUTDATED***
        """
        return False

    def rayIntersect(self, ray):
        """
        Returns the first intersection of a ray with the surface
        in parametric form (u,v,t) if such a solution exists.
        Otherwise, returns None.
        """
        return None

    def x(self, v, w):
        """
        Parametric equation for x
        """
        return 0 * u.cm

    def y(self, v, w):
        """
        Parametric equation for y
        """
        return 0 * u.cm

    def z(self, v, w):
        """
        Parametric equation for z
        """
        return 0 * u.cm

    def du(self, v, w):
        """
        First partial derivative with respect to u
        """
        return np.array((0, 0, 0))

    def dv(self, v, w):
        """
        First partial derivative with respect to v
        """
        return np.array((0, 0, 0))

    def getPoint(self, v, w):
        """
        Returns a point for parameters u and v
        """
        return u.Quantity((self.x(v, w), self.y(v, w), self.z(v, w)))

    def getNormal(self, v, w):
        """
        Returns the unit normal for parameters u and v
        """
        cross = np.cross(self.du(v, w), self.dv(v, w))
        return normalize(cross)
