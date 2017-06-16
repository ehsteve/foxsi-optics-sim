from foxsisim.module import Module
from foxsisim.detector import Detector
from foxsisim.source import Source
from foxsisim.plotting import scatterHist

import matplotlib.pyplot as plt

import astropy.units as u
import numpy as np

focalLength = 200.0
segmentLength = 30.0
radii = [5.15100, 4.90000, 4.65900, 4.42900, 4.21000, 4.00000, 3.79900]  # 7 shell radii
radii = [5.1510]
module = Module(seglen=segmentLength, focal=focalLength, radii=radii, conic=True, core_radius=0.1)

detector = Detector()
num_pixels = 256
pixel_size = 75 * u.micron
width = height = num_pixels * pixel_size
detector = Detector(width=width.to('cm').value,
                    height=height.to('cm').value,
                    normal=[0,0,1],
                    center=[0,0,200.0+30],
                    reso =[num_pixels, num_pixels])
source_distance = -10000000
offaxis_angle_arcmin = 0.0
center = [ source_distance * np.sin(np.deg2rad(offaxis_angle_arcmin/60.0)) , 0.0 , source_distance ]
print(center)
source = Source(type='atinf',
                center=center,
                color=[1,1,1])
#source = Source(type='atinf', center=center, color=[1,1,1])

number_of_rays = 100

#%timeit rays = source.generateRays(module.targetFront, number_of_rays)
rays = source.generateRays(module.targetFront, 10, grid=[10,10])
module.passRays(rays, robust=True)
detector.catchRays(rays)

one_bounce = [ray for ray in rays if ray.bounces == 2]
print(len(one_bounce))
module.plot2D(rays=rays)
plt.xlim(0, 200)
plt.ylim(3,6)
plt.show()
