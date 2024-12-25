from antenna import Antenna
from source import Source

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from typing import List
import imageio.v2 as imageio

EPS_0 = 8.854e-12
MU_0 = 4e-7 * np.pi

def max_except_source(E: np.ndarray, i: List[int], j: List[int]):
  mask = np.ones(E.shape, dtype=bool)
  for _i in i:
    for _j in j:
      mask[_i, _j] = False
  return np.max(E[mask])

def min_except_source(E: np.ndarray, i: List[int], j: List[int]):
  mask = np.ones(E.shape, dtype=bool)
  for _i in i:
    for _j in j:
      mask[_i, _j] = False
  return np.min(E[mask])
    
class ElectoMagneticMesh:
  def __init__(self, grid_size: int, dx: float, dy: float, antenna: Antenna = None, sources: List[Source] = []):
    self.grid_size = grid_size
    self.dx = dx
    self.dy = dy
    
    self.Ex = np.zeros((grid_size, grid_size))
    self.Ey = np.zeros((grid_size, grid_size))
    self.Hz = np.zeros((grid_size, grid_size))
    self.Ex_next = np.zeros((grid_size, grid_size))
    self.Ey_next = np.zeros((grid_size, grid_size))
    self.Hz_next = np.zeros((grid_size, grid_size))
    
    self.Ex_max = 0
    
    self.antenna = antenna
    self.sources = sources
    self.time = []
    self.max_field_coord = []
    self.max_field = []
    
  def add_antenna(self, antenna: Antenna):
    self.antenna = antenna
    
  def add_source(self, source: Source):
    self.sources.append(source)

  def check_collision(self, x: float, y: float):
    return y >= self.antenna.antenna_equation(x) > y - self.dy

  def check_inside_antenna(self, x: float, y: float):
    return y < self.antenna.antenna_equation(x)
  
  def renew_focus_field(self):
    focus_y = int(self.antenna.y0 / self.dy) + int(self.antenna.focal_length / self.dy)
    focus_x = int(self.antenna.x0 / self.dx) + self.grid_size//2
    self.Ex_max = max(self.Ex_max, (self.Ex[focus_y, focus_x])**2 + (self.Ey[focus_y, focus_x])**2)

  def renew_axe_field(self):
    I_axe = np.array((self.Ex[:, self.grid_size//2])**2 + (self.Ey[:, self.grid_size//2])**2)
    self.max_field_coord.append(np.argmax(I_axe) * self.dy)
    self.max_field.append(np.max(I_axe))

    # self.max_field_coord.append(np.argmax(abs(self.Ex[:self.grid_size//2, self.grid_size//2])) * self.dy)

    
  def _calculate_fields(self, dt: float = 0.0025 / (2 * 3e8)):
    for i in range(1, self.grid_size):
      for j in range(1, self.grid_size):
        self.Hz[i, j] -= (dt  / MU_0) * ((self.Ey[i, j] - self.Ey[i, j - 1]) / self.dx - (self.Ex[i, j] - self.Ex[i - 1, j]) / self.dy)

    for i in range(self.grid_size - 1):
      for j in range(self.grid_size - 1):
        if self.check_collision((j - self.grid_size//2) * self.dx, i * self.dy):
          self.Ex[i, j] = -1.0 * self.Ey[i, j] * self.antenna.tangent((j - self.grid_size//2) * self.dx)
          continue
        if self.check_inside_antenna((j - self.grid_size//2) * self.dx, i * self.dy):
          self.Ex[i, j] = 0
          self.Ey[i, j] = 0
          continue
        self.Ex[i, j] += (dt / EPS_0) * ((self.Hz[i + 1, j] - self.Hz[i, j]) / self.dy)
        self.Ey[i, j] -= (dt / EPS_0) * ((self.Hz[i, j + 1] - self.Hz[i, j]) / self.dx)

  def calculate(self, num_steps: int = 100, dt: float = 0.0025 / (2 * 3e8), visualise: bool = 'False', filename: str = "test"):
    frames = []
    for source in self.sources:
      self.Ex[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = source.source_func(0)
      self.Hz[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = np.sqrt(EPS_0 / MU_0) * source.source_func(0)

    for t in tqdm(range(num_steps)):
      self._calculate_fields()

      if t < num_steps // 5:
        for source in self.sources:
          self.Ex[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = source.source_func(t * dt)
          self.Hz[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = np.sqrt(EPS_0 / MU_0) * source.source_func(t * dt)
      else:
        #self.renew_focus_field()
        self.time.append(t * dt)
        self.renew_axe_field()

      if visualise:
        if t % 10 == 0:
          self.visualize(t * dt, frames)

    if visualise:
      imageio.mimsave(f'./visuals/{filename}.gif', frames, duration=0.4)


  def visualize(self, t: float, frames):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(0, self.grid_size * self.dx)
    ax.set_ylim(0, self.grid_size * self.dy)
    ax.set_xlabel("X (м)")
    ax.set_ylabel("Y (м)")
    ax.set_title("Электромагнитное поле")
    ax.axis("equal")

    if self.antenna:
      self.antenna.plot_antenna(ax=ax)

    idx = []
    jdx = []
    for source in self.sources:
      idx.append(int(source.y0 / self.dy))
      jdx.append(int(source.x0 / self.dx) + self.grid_size//2)
    min_E = min_except_source(self.Ex, idx, jdx)
    max_E = max_except_source(self.Ex, idx, jdx)
    edge = max(abs(min_E), abs(max_E))
    extent = [self.antenna.x0 - self.grid_size * self.dx / 2, self.antenna.x0 + self.grid_size * self.dx / 2, \
              self.antenna.y0, self.antenna.y0 + self.grid_size * self.dy]
    ax.imshow(self.Ex, cmap='hot', vmin = -edge, vmax = edge, extent=extent, origin='lower')
    ax.set_title(f"Время {round(t, 10)}")
    ax.set_xlabel("X (м)")
    ax.set_ylabel("Y (м)")
    ax.axis("equal")

    plt.savefig('temp.png')
    plt.close()
    frames.append(imageio.imread('temp.png'))

'''  def visualize(self, num_steps: int = 100, dt: float = 0.0025 / (2 * 3e8), filename: str = "test"):
    for source in self.sources: 
      self.Ex[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = source.source_func(0)
      self.Hz[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = np.sqrt(EPS_0 / MU_0) * source.source_func(0)

    for t in tqdm(range(num_steps)):
      self._calculate_fields()
          
      if t < num_steps // 5:
        for source in self.sources:
          self.Ex[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = source.source_func(t * dt)
          self.Hz[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = np.sqrt(EPS_0 / MU_0) * source.source_func(t * dt)
      else:
        self.renew_focus_field()

      if t % 10 == 0:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlim(0, self.grid_size * self.dx)
        ax.set_ylim(0, self.grid_size * self.dy)
        ax.set_xlabel("X (м)")
        ax.set_ylabel("Y (м)")
        ax.set_title("Электромагнитное поле")
        ax.axis("equal")

        if self.antenna:
          self.antenna.plot_antenna(ax=ax)

        idx = []
        jdx = []
        for source in self.sources:
          idx.append(int(source.y0 / self.dy))
          jdx.append(int(source.x0 / self.dx) + self.grid_size//2)
        min_E = min_except_source(self.Ex, idx, jdx)
        max_E = max_except_source(self.Ex, idx, jdx)
        edge = max(abs(min_E), abs(max_E))
        extent = [self.antenna.x0 - self.grid_size * self.dx / 2, self.antenna.x0 + self.grid_size * self.dx / 2, \
          self.antenna.y0, self.antenna.y0 + self.grid_size * self.dy]
        ax.imshow(self.Ex, cmap='hot', vmin = -edge, vmax = edge, extent=extent, origin='lower')
        ax.set_title(f"Время {round(t * dt, 10)}")
        ax.set_xlabel("X (м)")
        ax.set_ylabel("Y (м)")
        ax.axis("equal")
        
        plt.savefig('temp.png')
        plt.close()
        frames.append(imageio.imread('temp.png'))
        
    imageio.mimsave(f'./visuals/{filename}.gif', frames, duration=0.4)
'''

      