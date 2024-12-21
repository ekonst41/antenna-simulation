import numpy as np
import matplotlib.pyplot as plt
from antenna import Antenna
from tqdm import tqdm
import imageio.v2 as imageio

frames = []

epsilon_0 = 8.854e-12
mu_0 = 4e-7 * np.pi

class Source:
  def __init__(self, x0: float, y0: float, source_func):
    self.x0 = x0
    self.y0 = y0
    self.source_func = source_func
    
class ElectoMagneticMesh:
  def __init__(self, grid_size: int, dx: float, dy: float, antenna: Antenna = None, source: Source = None):
    self.grid_size = grid_size
    self.dx = dx
    self.dy = dy
    self.Ex = np.zeros((grid_size, grid_size))
    self.Ey = np.zeros((grid_size, grid_size))
    self.Hz = np.zeros((grid_size, grid_size))
    self.antenna = antenna
    self.source =source
    
  def add_antenna(self, antenna: Antenna):
    self.antenna = antenna
    
  def add_source(self, source: Source):
    self.source = source

  def check_collision(self, x: float, y: float):
    return y >= self.antenna.antenna_equation(x) > y - self.dy

  def check_inside_antenna(self, x: float, y: float):
    return y < self.antenna.antenna_equation(x)
    
  def _calculate_fields(self, dt: float = 0.0025 / (2 * 3e8)):
    for i in range(1, self.grid_size):
        for j in range(1, self.grid_size):
            self.Hz[i, j] -= (dt  / mu_0) * ((self.Ey[i, j] - self.Ey[i, j - 1]) / self.dx - (self.Ex[i, j] - self.Ex[i - 1, j]) / self.dy)

    for i in range(self.grid_size - 1):
        for j in range(self.grid_size - 1):
            if self.check_collision((j - self.grid_size//2) * self.dx, i * self.dy):
              self.Ex[i, j] = -1.0 * self.Ey[i, j] * self.antenna.tangent((j - self.grid_size//2) * self.dx)
              continue
            if self.check_inside_antenna((j - self.grid_size//2) * self.dx, i * self.dy):
              self.Ex[i, j] = 0
              self.Ey[i, j] = 0
              continue
            self.Ex[i, j] += (dt / epsilon_0) * ((self.Hz[i + 1, j] - self.Hz[i, j]) / self.dy)
            self.Ey[i, j] -= (dt / epsilon_0) * ((self.Hz[i, j + 1] - self.Hz[i, j]) / self.dx)
    # print(self.Hz, self.Ex, self.Ey)
            
  def visualize(self, num_steps: int = 100, dt: float = 0.0025 / (2 * 3e8)):
      self.Ex[int(self.source.y0 / self.dy), int(self.source.x0 / self.dx) + self.grid_size//2] = self.source.source_func(0)
      self.Hz[int(self.source.y0 / self.dy), int(self.source.x0 / self.dx) + self.grid_size//2] = np.sqrt(epsilon_0 / mu_0) * self.source.source_func(0)
      '''fig, ax = plt.subplots(figsize=(6, 6))
      ax.set_xlim(0, self.grid_size * self.dx)  # Установка масштаба по оси X
      ax.set_ylim(0, self.grid_size * self.dy)  # Установка масштаба по оси Y
      ax.set_xlabel("X (м)")
      ax.set_ylabel("Y (м)")
      ax.set_title("Электромагнитное поле")
      ax.axis("equal")'''

      for t in tqdm(range(num_steps)):
        self._calculate_fields()
          
        if t < num_steps:
          self.Ex[int(self.source.y0 / self.dy), int(self.source.x0 / self.dx) + self.grid_size//2] = self.source.source_func(t * dt)  # Применение источника
          self.Hz[int(self.source.y0 / self.dy), int(self.source.x0 / self.dx) + self.grid_size//2] = np.sqrt(epsilon_0 / mu_0) * self.source.source_func(t * dt)
          #print(self.Ex[int(self.source.y0 / self.dy), int(self.source.x0 / self.dx) + self.grid_size//2])
        '''if t == num_steps // 10:
          self.Ex[int(self.source.y0 / self.dy), int(self.source.x0 / self.dx) + self.grid_size//2] = 0
        '''
        if t % 10 == 0:
          #print(np.min(self.Ex), np.max(self.Ex))
          #print(self.Ex[0, 0])
          fig, ax = plt.subplots(figsize=(6, 6))
          ax.set_xlim(0, self.grid_size * self.dx)  # Установка масштаба по оси X
          ax.set_ylim(0, self.grid_size * self.dy)  # Установка масштаба по оси Y
          ax.set_xlabel("X (м)")
          ax.set_ylabel("Y (м)")
          ax.set_title("Электромагнитное поле")
          ax.axis("equal")
          #ax.clear()

          if self.antenna:
            self.antenna.plot_antenna(ax=ax)

                # Отрисовка электрического поля
          ax.imshow(self.Ex, cmap='hot', extent= [self.antenna.x0 - self.grid_size * self.dx / 2, self.antenna.x0 + self.grid_size * self.dx / 2, self.antenna.y0, self.antenna.y0 + self.grid_size * self.dy], origin='lower')
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
      imageio.mimsave('./visuals/antenna3.gif', frames, duration=0.4)  # duration - время отображения каждого кадра в секундах

#plt.show()
      
      