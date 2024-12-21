from antenna import Antenna
from source import Source

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from typing import List
import imageio.v2 as imageio

frames = []

epsilon_0 = 8.854e-12
mu_0 = 4e-7 * np.pi

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
    
    self.antenna = antenna
    self.sources = sources
    
  def add_antenna(self, antenna: Antenna):
    self.antenna = antenna
    
  def add_source(self, source: Source):
    self.sources.append(source)

  def check_collision(self, x: float, y: float):
    return y >= self.antenna.antenna_equation(x) > y - self.dy

  def check_inside_antenna(self, x: float, y: float):
    return y < self.antenna.antenna_equation(x)
    
  def _calculate_fields(self, dt: float = 0.0025 / (2 * 3e8)):
    for i in range(1, self.grid_size):
        for j in range(1, self.grid_size):
            self.Hz_next[i, j] -= (dt  / mu_0) * ((self.Ey[i, j] - self.Ey[i, j - 1]) / self.dx - (self.Ex[i, j] - self.Ex[i - 1, j]) / self.dy)

    for i in range(self.grid_size - 1):
      for j in range(self.grid_size - 1):
        if self.check_collision((j - self.grid_size//2) * self.dx, i * self.dy):
          self.Ex_next[i, j] = -1.0 * self.Ey[i, j] * self.antenna.tangent((j - self.grid_size//2) * self.dx)
          continue
        if self.check_inside_antenna((j - self.grid_size//2) * self.dx, i * self.dy):
          self.Ex_next[i, j] = 0
          self.Ey_next[i, j] = 0
          continue
        self.Ex_next[i, j] += (dt / epsilon_0) * ((self.Hz[i + 1, j] - self.Hz[i, j]) / self.dy)
        self.Ey_next[i, j] -= (dt / epsilon_0) * ((self.Hz[i, j + 1] - self.Hz[i, j]) / self.dx)
        
        self.Ex = np.copy(self.Ex_next)
        self.Ey = np.copy(self.Ey_next)
        self.Hz = np.copy(self.Hz_next)
            
  def visualize(self, num_steps: int = 100, dt: float = 0.0025 / (2 * 3e8)):
    for source in self.sources: 
      self.Ex[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = source.source_func(0)
      self.Hz[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = np.sqrt(epsilon_0 / mu_0) * source.source_func(0)

    for t in tqdm(range(num_steps)):
      self._calculate_fields()
          
      if t < num_steps // 5:
        for source in self.sources:
          self.Ex[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = source.source_func(t * dt)  # Применение источника
          self.Hz[int(source.y0 / self.dy), int(source.x0 / self.dx) + self.grid_size//2] = np.sqrt(epsilon_0 / mu_0) * source.source_func(t * dt)

      if t % 10 == 0:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlim(0, self.grid_size * self.dx)  # Установка масштаба по оси X
        ax.set_ylim(0, self.grid_size * self.dy)  # Установка масштаба по оси Y
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
          # extent= [self.antenna.x0 - self.grid_size * self.dx / 2, self.antenna.x0 + self.grid_size * self.dx / 2, self.antenna.y0, self.antenna.y0 + self.grid_size * self.dy]
        ax.set_title(f"Время {round(t * dt, 10)}")
        ax.set_xlabel("X (м)")
        ax.set_ylabel("Y (м)")
        ax.axis("equal")

          # Обновление графика
          #plt.draw()
          #plt.pause(0.01)  # Пауза для обновления
        plt.savefig('temp.png')
        plt.close()
        frames.append(imageio.imread('temp.png'))
        
    imageio.mimsave('./visuals/antenna4.gif', frames, duration=0.4)  # duration - время отображения каждого кадра в секундах

#plt.show()
      
      