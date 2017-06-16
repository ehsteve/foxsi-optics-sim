"""
Created on Jul 8, 2011

@author: rtaylor
"""
import numpy as np
import astropy.units as u
from foxsisim.mymath import to_cm


class Ray:
    """
    Holds the position and orientation vectors of a light ray
    """
    @u.quantity_input(pos=u.cm, ori=u.cm, src=u.cm, des=u.cm)
    def __init__(self,
                 pos=[0, 0, 0]*u.cm,
                 ori=[0, 0, 1]*u.cm,
                 src=[0, 0, 0]*u.cm,
                 des=[0, 0, 0]*u.cm,
                 dead=False,
                 tag=None,
                 bounces=0,
                 energy=None):
        """
        Constructor

        Parameters:
            pos:     current position coordinate
            ori:     current orientation vector, normalized
            src:     original source coordinate
            des:     final destination coordinate
            dead:    is the ray no longer moving
            tag:     used to identify the ray in some way, can store pointer to
                     source object
            bounces: the number of times the ray has hit a surface
            energy:  the energy of the ray (in keV)
            num:     a number to identify the ray
            hist:    the history of the ray, a list of positions of all the interactions
        """
        self.pos = pos
        self.ori = ori
        self.src = src
        self.des = des
        self.dead = dead
        self.tag = tag
        self.bounces = bounces
        self.energy = energy
        self.num = 0
        self.hist = []

    def inRange(self, t):
        """
        Is the t parameter in the desired range?
        """
        if t > 0.0000000000001:
            return True
        else:
            return False

    def x(self, t):
        """
        Parametric equation for x
        """
        return self.pos[0] + self.ori[0] * t

    def y(self, t):
        """
        Parametric equation for y
        """
        return self.pos[1] + self.ori[1] * t

    def z(self, t):
        """
        Parametric equation for z
        """
        return self.pos[2] + self.ori[2] * t

    def getPoint(self, t):
        """
        Returns a point for parameter t
        """
        return u.Quantity([self.x(t), self.y(t), self.z(t)])

    def plot3D(self, axes, color='r'):
        """
        Generates a 3d plot of the ray in the given figure
        """
        p1 = self.getPoint(0)
        p2 = self.getPoint(1)
        axes.plot3D(to_cm((p1[0], p2[0])), to_cm((p1[1], p2[1])),
                    to_cm((p1[2], p2[2])), '-' + color)

    def moveToZ(self, z):
        """
        Move ray along ori vector to a specific z value.
        Returns true on success, else returns false.
        """
        # if already there
        if self.pos[2] == z:
            return True
        # if ray cannot move in z direction
        if self.ori[2] == 0:
            return False
        # find t
        t = (z - self.pos[2]) / self.ori[2]
        # if desired z is in opposite direction
        if t < 0:
            return False
        # move ray
        self.pos = self.getPoint(t)
        return True
