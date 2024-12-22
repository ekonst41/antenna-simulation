from antenna import Antenna
from mesh import ElectoMagneticMesh, Source

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

grid_size = 200
dx, dy = 0.0025, 0.0025

def source_function(t):
  return 0.1 * np.cos(2 * np.pi * 3.0e8 * t)


def field_in_focus(min: float = 0.05, max: float = 1, num: int = 5):
  fields = {}
  focuses = np.linspace(min, max, num)
  
  sources = []
  for i in range(grid_size):
    source = Source((i - grid_size // 2) * dx , 0.2, source_function)
    sources.append(source)
  
  for focus_lenght in tqdm(focuses):
    antenna = Antenna(0, 0, focus_lenght, focus_lenght)

    mesh = ElectoMagneticMesh(grid_size, dx, dy, sources=sources)
    mesh.add_antenna(antenna)
    #mesh.visualize(500)
    mesh.calculate(num_steps=500, visualise=True, filename='value_in_focus')
    fields[focus_lenght] = mesh.Ex_max
  return fields

fields = field_in_focus(max=0.4, num=15)


focus_lengths = list(fields.keys())
field_values = list(fields.values())

df = pd.DataFrame({
  'focus_lengths': focus_lengths,
  'field_values': field_values,
})

df.to_csv('data/field_in_focus_005-040.csv', index=False)


'''# Построение графика
plt.figure(figsize=(10, 6))
plt.plot(focus_lengths, field_values, marker='o', linestyle='-', color='b')
plt.title('Зависимость поля в фокусе от фокусного расстояния')
plt.xlabel(r'$ f $')
plt.ylabel(r'$ E $')
plt.grid(True)
plt.show()'''


