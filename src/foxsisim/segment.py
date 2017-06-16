"""
Created on Jul 8, 2011

@author: rtaylor
"""
from foxsisim.surface import Surface
from foxsisim.mymath import to_cm, normalize
import numpy as np
from scipy.optimize import fsolve
import astropy.units as u


class Segment(Surface):
    """
    A shell segment surface defined by parametric equations. It is aligned
    such that its axis is in the z direction, its wide end is facing down
    (negative z) and its small end is facing up (positive z). Usually, the
    radius of the small end (r1) is computed by the constructor. However,
    the user may supply a value for r1, in which case the wide radius (r0)
    is automatically recomputed. When changing segment dimensions after
    instantiation, use the updateDims method.
    """

    @u.quantity_input(base=u.cm, seglen=u.cm, ang=u.deg, r0=u.cm)
    def __init__(self,
                 base=[0, 0, 0]*u.cm,
                 seglen=30.0 * u.cm,
                 ang=0.006 * u.deg,
                 r0=5.5 * u.cm,
                 r1=None
                 ):
        """
        Constructor

        Parameters:
            base:    the center point of the wide end of the segment
            seglen:  the axial length of the segment
            ang:     angle of the segment side to the axis
            r0:      radius of the wide end of the segment
            r1:      radius of the small end of the segment
        """
        # instantiate
        Surface.__init__(self)
        self.base = base
        self.seglen = seglen
        self.ang = ang
        self.updateDims(r0, r1)

    @u.quantity_input(r0=u.cm)
    def updateDims(self, r0, r1=None):
        """
        Takes an r0 value and sets the appropriate r1 and ch values (ch is an
        important internal parameter). Alternatively, we can pass r1 and it
        sets the appropriate r0 and ch.
        """
        # update r0 and r1
        if r1 is None:
            self.r0 = r0
            self.r1 = self.r0 - self.seglen * np.tan(self.ang)
        else:
            self.r1 = r1
            self.r0 = self.r1 + self.seglen * np.tan(self.ang)
        # calculate cone height
        self.ch = self.seglen / (1 - self.r1 / self.r0) # cone height is used by the parametric eqns
        # check vars
        if self.r0 <= 0 or self.r0 < self.r1 or self.r1 <= 0:
            print('error: invalid segment dimensions')

    def inRange(self, u, v):
        """
        Are the u and v parameters in the desired range?
        """
        if u >= 0 and u <= self.seglen:
            return True
        else:
            return False

    def existsInOctant(self,octant):
        """
        Returns whether the surface exists in the supplied octant ***OUTDATED***
        """

        def system(params):
            """
            Returns the displacement of a given point on the surface with
            the nearest point in the octant. Requires 3 parameters, but
            only two are used: (u,v,_)
            """
            pnt1 = self.getPoint(params[0], params[1])
            pnt2 = octant.nearestPoint(pnt1)
            return (pnt1[0]-pnt2[0], pnt1[1]-pnt2[1], pnt1[2]-pnt2[2])

        limits = octant.getLimits()
        ctrx = (limits[0][0]+limits[0][1])/2
        ctry = (limits[1][0]+limits[1][1])/2
        guess = (0, np.arctan2(ctry, ctrx),0)
        X, infodict, ier, mesg = fsolve(system, guess, full_output=1) #@UnusedVariable
        valid = False

        debug = False
        # if solution was found
        if(ier == 1):
            # fsolve sometimes returns nonexistant solutions when the octant is in the interior of the surface
            Y = system(X)
            if abs(Y[0])+abs(Y[1])+abs(Y[2]) > 0.0001:
                if debug: print('Solution found, but is nonzero')
            elif not self.inRange(X[0], X[1]):
                if debug: print('Solution found, but outside of surface\'s parameter range:')
            else:
                valid = True
                if debug:
                    print('Valid solution found:')
            if debug:
                print('[u,v,_] = '), X
                print('[x,y,z] = '), self.getPoint(X[0],X[1])
        # no solution found
        else:
            if debug:
                print('No solution found: '), mesg
            if debug:
                print('Last iteration: '), X

        return valid

    def rayIntersect(self, ray):
        """
        Returns the first intersection of a ray with the surface
        in parametric form (u,v,t) if such a solution exists.
        Otherwise, returns None.
        """
        # display debug messages
        debug = False

        def system(params):
            """
            Returns the displacement from a point on the ray to a
            point on the surface. Requires 3 parameters: (u,v,t)
            """
            x = self.x(params[0],params[1]) - ray.x(params[2])
            y = self.y(params[0],params[1]) - ray.y(params[2])
            z = self.z(params[0],params[1]) - ray.z(params[2])
            return (x, y, z)

        # make sure ori is a unit vector
        ray.ori = normalize(ray.ori)

        # solve with guess at current ray position
        guess = u.Quantity(ray.pos[2] - self.base[2], np.arctan2(ray.pos[1], ray.pos[0]), 0).value
        X, infodict, ier, mesg = fsolve(system, guess, full_output=1) #@UnusedVariable

        # if solution found
        if ier == 1:
            # if solution is out of range
            if not ray.inRange(X[2]):

                # debug message
                if debug: print('Segment: first solution out of range. Trying again.')

                # ray can intersect at most twice, so create a second guess
                tsup = 2*self.seglen
                pos = ray.getPoint(tsup) # a point really far away from current position
                guess = (pos[2]-self.base[2],atan2(pos[1],pos[0]),tsup)
                X,infodict,ier,mesg = fsolve(system,guess,full_output=1) #@UnusedVariable

        # return solution if exists
        result = None
        if(ier == 1):
            if not self.inRange(X[0],X[1]):
                if debug: print('Segment: final solution out of surface parameter range:')
            elif not ray.inRange(X[2]):
                if debug: print('Segment: final solution out of ray parameter range:')
            else:
                if debug: print('Segment: valid solution found:')
                result = X
            if debug:
                print('   [u,v,t] = '),X
                print('   [x,y,z] = '),self.getPoint(X[0],X[1])
        # no solution found
        else:
            if debug:
                print('Segment: no solution found: '),mesg
                print('   [u,v,t] = '),X

        return result

    def x(self, v, w):
        """
        Parametric equation for x
        """
        return self.r0 * np.cos(w) * (self.ch - v) / self.ch + self.base[0]

    def y(self, v, w):
        """
        Parametric equation for y
        """
        return self.r0 * np.sin(w) * (self.ch - v) / self.ch + self.base[1]

    def z(self, v, w):
        """
        Parametric equation for z
        """
        return v + self.base[2]

    def du(self, v, w):
        """
        First partial derivative with respect to u
        """
        dx = -self.r0 * np.cos(w) / self.ch
        dy = -self.r0 * np.sin(w) / self.ch
        dz = 1
        return np.array((dx,dy,dz))

    def dv(self, v, w):
        """
        First partial derivative with respect to v
        """
        dx = -self.r0 * np.sin(w) * (self.ch - v) / self.ch
        dy = self.r0 * np.cos(w) * (self.ch - v) / self.ch
        dz = 0
        return np.array((dx,dy,dz))

    def plot2D(self, axes, color = 'b'):
        """
        Plots a 2d cross section of the segment. All units in cm.
        """
        x = u.Quantity((self.base[2], (self.base[2]+self.seglen))).to('cm').value
        y = u.Quantity((self.r0, self.r1)).to('cm').value
        axes.plot(x, y, '-' + color)
        y = u.Quantity((-self.r0, -self.r1)).to('cm').value
        axes.plot(x, y, '-' + color)

    def plot3D(self, axes, color = 'b'):
        """
        Generates a 3d plot of the segment in the given figure
        """
        v = np.linspace(0, 2 * np.pi,20)
        for i in range(len(v)):
            p1 = self.getPoint(0, v[i])
            p2 = self.getPoint(0, v[i-1])
            p3 = self.getPoint(self.seglen, v[i-1])
            p4 = self.getPoint(self.seglen, v[i])
            axes.plot3D(to_cm([p1[0], p2[0], p3[0], p4[0]]),
                        to_cm([p1[1], p2[1], p3[1], p4[1]]),
                        to_cm([p1[2], p2[2], p3[2], p4[2]]),
                        '-'+color)

    def targetFront(self,a,b):
        """
        Takes two list arguments of equal size, the elements of which range from 0 to 1.
        Returns an array of points that exist on the circle defined by the wide end of
        the segment.
        """
        n = len(a)
        pnts = np.zeros((n,3))
        for i in range(n):
            ang = 2 * 180*u.deg * b[i]
            r = self.r0 * np.sqrt(a[i]) # we want the polar coords to have a uniform cartesian distribution
            pnts[i,0] = self.base[0] + r * np.cos(ang)
            pnts[i,1] = self.base[1] + r * np.sin(ang)
            pnts[i,2] = self.base[2]
        return pnts

    def targetBack(self,a,b):
        """
        Takes two list arguments of equal size, the elements of which range from 0 to 1.
        Returns an array of points that exist on the circle defined by the small end of
        the segment.
        """
        n = len(a)
        pnts = np.zeros((n,3))
        for i in range(n):
            ang = 2 * 180 * u.deg * b[i]
            r = self.r1 * np.sqrt(a[i]) # we want the polar coords to have a uniform cartesian distribution
            pnts[i][0] = self.base[0] + r * np.cos(ang)
            pnts[i][1] = self.base[1] + r * np.sin(ang)
            pnts[i][2] = self.base[2] + self.seglen
        return pnts
