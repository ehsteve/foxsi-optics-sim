"""
Created on Jul 19, 2011

@author: rtaylor
"""
from foxsisim.shell import Shell
from foxsisim.circle import Circle
from foxsisim.mymath import reflect, calcShellAngle, to_cm, normalize, norm
import astropy.units as u
import numpy as np


class Module:
    """
    A complete foxsi module. By default, it consists of seven nested shells.
    """
    @u.quantity_input(base=u.cm, seglen=u.cm, focal=u.cm, radii=u.cm)
    def __init__(self,
                 base=[0, 0, 0]*u.cm,
                 seglen=30.0*u.cm,
                 focal=200.0*u.cm,
                 radii=[5.151, 4.9, 4.659, 4.429, 4.21, 4.0, 3.799]*u.cm,
                 conic=True,
                 shield=True,
                 core_radius=None
                 ):
        """
        Constructor

        Parameters:
            base:       the center point of the wide end of the segment
            seglen:     the axial length of each segment
            focal:      the focal length, measured from the center of the module
            radii:      a list of radii, one for each shell from biggest to
                        smallest
            angles:     optional parameter to overwrite the shell angles computed
                        by constructor
            conic:      if True, use a conic approximation to the Wolter-I parabola-hyperbola.
            shield:     if True, create a shield at the core of the optic to block straight-thru
                        flux.
            shield_radius If shield is True then use this value for the shield radius
        """
        angles = calcShellAngle(radii, focal)

        self.shells = []
        for i, r in enumerate(radii):
            self.shells.append(Shell(base=base, focal=focal, seglen=seglen, ang=angles[i], r=r, conic=conic))

        # inner core (blocks rays going through center of module)
        # not sure if the core must have an angle to it so made it very small
        # and made coreFaces match in radius
        if shield is True:
            self.shield = True
            if core_radius is None:
                r0 = self.shells[-1].front.r0
                r1 = r0 - seglen * np.tan(4 * angles[-1])
                #ang = atan((r0 - r1) / (2 * seglen))
            else:
                r0 = core_radius
                r1 = core_radius
            #self.core = Segment(base=base, seglen=2 * seglen, ang=0, r0=r0)
            self.coreFaces = [Circle(center=base, normal=[0, 0, 1], radius=r0),
                              Circle(center=u.Quantity([base[0], base[1],
                                             base[2] + 2 * seglen]),
                                     normal=[0, 0, -1], radius=r0)]
        else:
            self.shield = None

    def getDims(self):
        """
        Returns the module's dimensions:
        [radius at wide end, radius at small end, length]
        """
        front = self.shells[0].front
        back = self.shells[0].back
        return [front.r0, back.r1, front.seglen + back.seglen]

    def getSurfaces(self):
        """
        Returns a list of surfaces
        """
        surfaces = []
        for shell in self.shells:
            surfaces.extend(shell.getSurfaces())
        #surfaces.append(self.core)
        surfaces.extend(self.coreFaces)
        return(surfaces)

    def passRays(self, rays, robust=False):
        """
        Takes an array of rays and passes them through the front end of
        the module.
        """
        # print('Module: passing ',len(rays),' rays')

        # get all module surfaces
        allSurfaces = self.getSurfaces()
        allSurfaces.remove(self.coreFaces[0])  # we'll test these seperately
        allSurfaces.remove(self.coreFaces[1])

        # create regions consisting of adjacent shells
        regions = [None for shell in self.shells]
        for i, shell in enumerate(self.shells):
            # innermost shell
            if i == len(self.shells) - 1:
                regions[i] = shell.getSurfaces()
                #regions[i].append(self.core)
            else:
                # outer shell (reflective side facing region)
                regions[i] = shell.getSurfaces()
                # nested shell (non reflective)
                regions[i].extend(self.shells[i + 1].getSurfaces())

        for ray in rays:
            if self.shield is not None:
                # skip rays that hit a core face
                if ray.pos[2] < self.coreFaces[0].center[2]:
                    sol = self.coreFaces[0].rayIntersect(ray)
                    if sol is not None:
                        sol = sol * u.cm
                        #print("ray hit face 0")
                        ray.pos = ray.getPoint(sol[2])
                        ray.bounces += 1
                        ray.dead = True
                        ray.des = ray.pos
                        ray.hist.append(ray.pos)
                        continue
                    else:
                        ray.moveToZ(self.coreFaces[0].center[2])
                elif ray.pos[2] > self.coreFaces[1].center[2]:
                    sol = self.coreFaces[1].rayIntersect(ray)
                    if sol is not None:
                        #print("ray hit face 1")
                        ray.pos = ray.getPoint(sol[2])
                        ray.bounces += 1
                        ray.dead = True
                        ray.des = ray.pos
                        ray.hist.append(ray.pos)
                        continue
                    else:
                        ray.moveToZ(self.coreFaces[1].center[2])

            # reset surfaces
            surfaces = [s for s in allSurfaces]
            firstBounce = True  # used for optimization

            # while ray is inside module
            while True:
                # find nearest ray intersection
                bestDist = None
                bestSol = None
                bestSurf = None
                for surface in surfaces:

                    sol = surface.rayIntersect(ray)
                    if sol is not None:
                        dist = norm(ray.getPoint(sol[2]) - ray.pos)
                        if bestDist is None or dist < bestDist:
                            bestDist = dist
                            bestSol = sol
                            bestSurf = surface

                # if a closest intersection was found
                if bestSol is not None:

                    # update ray
                    ray.pos = ray.getPoint(bestSol[2])
                    ray.hist.append(ray.pos)
                    ray.bounces += 1
                    #print("%i ray bounce number %i" % (ray.num, ray.bounces))

                    x = reflect(ray.ori,
                                bestSurf.getNormal(bestSol[0], bestSol[1]),
                                ray.energy)
                    #print(x)
                    # if reflected
                    if x is not None:
                        # update ori to unit vector reflection
                        ray.ori = normalize(x)
                    # otherwise, no reflection means ray is dead
                    else:
                        ray.dead = True
                        #print("%i ray killed by reflect" % ray.num)
                        break

                    # knowing the surface it has just hit, we can
                    # narrow down the number of surface to test

                    # remove shells the ray cannot even 'see'
                    if firstBounce:
                        firstBounce = False
                        for region in regions:
                            if bestSurf is region[0] or bestSurf is region[1]:
                                surfaces = [s for s in region]
                                break

                    # assuming each segment can be hit no more than once
                    # eliminate the surface from our list
                    if not robust:
                        surfaces.remove(bestSurf)

                # if no intersection, ray can exit module
                else:
                    break

    def plot2D(self, axes, color='b'):
        """
        Plots a 2d cross section of the module
        """
        for shell in self.shells:
            shell.plot2D(axes, color)
        if self.shield is not None:
            # plot core
            #self.core.plot2D(axes, color)
            #base = self.core.base
            r0 = self.coreFaces[0].radius
            r1 = self.coreFaces[1].radius
            z0 = self.coreFaces[0].center[2]
            z1 = self.coreFaces[1].center[2]
            #seglen = self.core.seglen

            axes.plot(to_cm((z0, z0)), to_cm((r0, -r0)), '-' + color)
            axes.plot(to_cm((z1, z1)), to_cm((r1, -r1)), '-' + color)

    def plot3D(self, axes, color='b'):
        """
        Generates a 3d plot of the module in the given figure
        """
        for shell in self.shells:
            shell.plot3D(axes, color)

    def targetFront(self, a, b):
        """
        Takes two list arguments of equal size, the elements of which range
        from 0 to 1. Returns an array of points that exist on the circle
        defined by the wide end of the module.
        """
        # must modify 'a' so that we dont return points from the core
        r0 = self.shells[0].front.r0
        r1 = self.coreFaces[0].radius
        a0 = (r1 / r0) ** 2  # the 'a' value that gives r1=sqrt(a)*r0
        adiff = 1 - a0
        for i in range(len(a)):
            a[i] = a[i] * adiff + a0
        return self.shells[0].targetFront(a, b)

    def targetBack(self, a, b):
        """
        Takes two list arguments of equal size, the elements of which range
        from 0 to 1. Returns an array of points that exist on the circle
        defined by the small end of the module.
        """
        # must modify 'a' so that we dont return points from the core
        r0 = self.shells[0].back.r1
        r1 = self.coreFaces[0].radius
        a0 = (r1 / r0) ** 2  # the 'a' value that gives r1=sqrt(a)*r0
        adiff = 1 - a0
        for i in range(len(a)):
            a[i] = a[i] * adiff + a0
        return self.shells[0].targetBack(a, b)
