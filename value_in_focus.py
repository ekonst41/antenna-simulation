from antenna import Antenna
from mesh import ElectoMagneticMesh, Source

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

grid_size = 200
dx, dy = 0.005, 0.005

def source_function(t):
  return 0.1 * np.cos(2 * np.pi * 3.0e8 * t)

def field_in_axe(focal_length: float):
  sources = []
  for i in range(1, grid_size):
    source = Source((i - grid_size // 2) * dx , 0.24, source_function)
    sources.append(source)
  antenna = Antenna(x0=0, y0=0, diameter=0.4, focal_length=focal_length)

  mesh = ElectoMagneticMesh(grid_size, dx, dy, sources=sources)
  mesh.add_antenna(antenna)
  mesh.calculate(num_steps=700, visualise=False)
  time1 = mesh.time
  max_field_coord1 = mesh.max_field_coord
  max_field1 = mesh.max_field
  return time1, max_field_coord1, max_field1


def field_in_focus(min: float = 0.05, max: float = 1, num: int = 5):
  fields = {}
  focuses = np.linspace(min, max, num)
  
  sources = []
  for i in range(grid_size):
    source = Source((i - grid_size // 2) * dx , 0.24, source_function)
    sources.append(source)
  
  for focus_lenght in tqdm(focuses):
    antenna = Antenna(x0=0, y0=0, diameter=0.4, focal_length=focus_lenght)

    mesh = ElectoMagneticMesh(grid_size, dx, dy, sources=sources)
    mesh.add_antenna(antenna)
    #mesh.visualize(500)
    mesh.calculate(num_steps=700, visualise=False)
    fields[focus_lenght] = mesh.Ex_max
  return fields

'''fields = field_in_focus(max=0.4, num=15)


focus_lengths = list(fields.keys())
field_values = list(fields.values())

df = pd.DataFrame({
  'focus_lengths': focus_lengths,
  'field_values': field_values,
})
'''
min = 0.1
max = 0.2
num = 21
focuses = np.linspace(min, max, num)
dt = 0.0025 / (2 * 3e8)
index = int(1e-9 / dt)

max_field = []
max_field_coord = []

for focus in tqdm(focuses):
  time, coord, field = field_in_axe(focus)
  max_field.append(np.max(field))
  max_field_coord.append(coord[np.argmax(field)])


'''df = pd.DataFrame({
  'time': time,
  'coord': coord,
  'field': field,
})
df.to_csv('data/field_on_axe_040_focus_length.csv', index=False)
'''
#max_field.append(np.max(field))
#max_field_coord.append(coord[np.argmax(field)])

df = pd.DataFrame({
    'focus': focuses,
    'real_coord': max_field_coord,
    'max_field': max_field,
  })
df.to_csv('data/field_on_axe_different_focals_01-02.csv', index=False)