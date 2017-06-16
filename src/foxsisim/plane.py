"""
Created on Jul 27, 2011

@author: rtaylor
"""
import numpy as np
import astropy.units as u
from foxsisim.surface import Surface
from foxsisim.mymath import normalize, to_cm, norm


class Plane(Surface):
    """
    A parallelogram surface defined by two axis vectors and an origin point.
    """

    @u.quantity_input(origin=u.cm, ax1=u.cm, ax2=u.cm)
    def __init__(self, origin=[-1, -1, 0] * u.cm, ax1=[2, 0, 0] * u.cm, ax2=[0, 2, 0] * u.cm):
        """
        Constructor

        Parameters:
            ax1:     first edge of parallelogram as a vector
            ax2:     second edge of parallelogram as a vector
            origin:  the origin coordinate for both axes
        """
        Surface.__init__(self)
        self.origin = origin
        self.ax1 = ax1
        self.ax2 = ax2

    def x(self, v, w):
        """
        Parametric equation for x
        """
        return self.origin[0] + v * self.ax1[0] + w * self.ax2[0]

    def y(self, v, w):
        """
        Parametric equation for y
        """
        return self.origin[1] + v * self.ax1[1] + w * self.ax2[1]

    def z(self, v, w):
        """
        Parametric equation for z
        """
        return self.origin[2] + v * self.ax1[2] + w * self.ax2[2]

    def getPoint(self, v, w):
        """
        Returns a point for parameters u and v
        """
        return self.origin + v * self.ax1 + w * self.ax2

    def du(self, u, v):
        """
        First partial derivative with respect to u
        """
        return self.ax1

    def dv(self, v, w):
        """
        First partial derivative with respect to v
        """
        return self.ax2

    def getWidth(self):
        """
        Returns the length of the parallelogram's first axis
        """
        return norm(self.ax1)

    def getHeight(self):
        """
        Returns the length of the parallelogram's second axis
        """
        return norm(self.ax2)

    def inRange(self, v, w):
        """
        Are the u and v parameters in the desired range?
        """
        if v >= 0 and v <= 1 and w >= 0 and w <= 1:
            return True
        else: return False

    def rayIntersect(self, ray):
        """
        Returns the first intersection of a ray with the surface
        in parametric form (u,v,t) if such a solution exists.
        Otherwise, returns None.
        """
        a = np.transpose(u.Quantity((self.ax1, self.ax2, -ray.ori)))
        b = np.transpose(ray.pos - self.origin)
        try: x = np.linalg.solve(a, b)
        except np.linalg.LinAlgError:
            # no solution exists (ray is parallel to plane)
            x = None
        if x is not None:
            if not self.inRange(x[0], x[1]):
                x = None
        return x

    def grid(self, a, b):
        """
        Takes two list arguments of equal size, the elements of which range
        from 0 to 1. Returns an array of points that exist at the corresponding
        locations on the parallelogram.
        """
        n = len(a)
        pnts = np.zeros((n,3)) * u.cm
        for i in range(n): pnts[i,:] = self.getPoint(a[i], b[i])
        return pnts

    def plot3D(self, axes, color = 'b'):
        """
        Generates a 3d plot of the plane in the given figure
        """
        p1 = self.origin
        p2 = self.getPoint(1,0)
        p3 = self.getPoint(1,1)
        p4 = self.getPoint(0,1)
        axes.plot3D(to_cm([p1[0],p2[0],p3[0],p4[0],p1[0]]),
                    to_cm([p1[1],p2[1],p3[1],p4[1],p1[1]]),
                    to_cm([p1[2],p2[2],p3[2],p4[2],p1[2]]),
                    '-'+color)
