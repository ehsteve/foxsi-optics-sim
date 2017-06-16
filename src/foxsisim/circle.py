"""
Created on Jul 27, 2011

@author: rtaylor
"""
from foxsisim.plane import Plane
import numpy as np
from numpy.linalg import norm
from foxsisim.mymath import normalize
import astropy.units as u


class Circle(Plane):
    """
    A circular surface defined by an origin point, a normal, and a radius.
    """

    @u.quantity_input(center=u.cm, radius=u.cm)
    def __init__(self, center=[0, 0, 0]*u.cm, normal=[0, 0, 1], radius=1*u.cm):
        """
        Constructor

        Parameters:
            center: center location of circle
            normal: surface normal of circle
            radius: radius of circle
        """
        # normal should be length 1
        normal = normalize(normal)

        # create rectangular dimensions
        if normal[0] == 0 and normal[2] == 0:  # normal is in y direction
            sign = normal[1]  # 1 or -1
            ax1 = radius * sign * np.array((0, 0, 1))
            ax2 = radius * np.array((0, 1, 0))
        else:
            ax1 = radius * np.cross(np.array([0, 1, 0]), normal)  # parallel to xz-plane
            ax2 = radius * np.cross(normal, ax1)

        Plane.__init__(self, origin=center, ax1=ax1, ax2=ax2)
        self.center = center
        self.normal = normal
        self.radius = radius

    def inRange(self, v, w):
        """
        Are the u and v parameters in the desired range?
        """
        if norm(self.getPoint(v, w).to('cm') - self.center.to('cm')) * u.cm <= self.radius:
            return True
        else:
            return False

    def plot3D(self, axes, color='b'):
        """
        Generates a 3d plot of the plane in the given figure
        """
        pass
