from antenna import Antenna
from mesh import ElectoMagneticMesh, Source

import numpy as np

grid_size = 200
dx, dy = 0.0025, 0.0025

def source_function(t):
    return 0.1 * np.cos(2 * np.pi * 3.0e8 * t)

antenna = Antenna(0, 0, 0.3, 0.25)

mesh = ElectoMagneticMesh(grid_size, dx, dy)

for i in range(grid_size):
  source = Source((i - grid_size // 2) * dx , 0.25, source_function)
  mesh.add_source(source)
  
mesh.add_antenna(antenna)


mesh.visualize(500)
