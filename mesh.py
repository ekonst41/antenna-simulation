import numpy as np
import matplotlib.pyplot as plt
from antenna import Antenna
from tqdm import tqdm

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
    
  def _calculate_fields(self, dt: float = 0.5):
    for i in range(1, self.grid_size):
        for j in range(1, self.grid_size):
            self.Hz[i, j] += (dt / (mu_0 * self.dx)) * (self.Ey[i, j] - self.Ey[i-1, j]) - \
                        (dt / (mu_0 * self.dy)) * (self.Ex[i, j] - self.Ex[i, j-1])

    for i in range(self.grid_size - 1):
        for j in range(self.grid_size - 1):
            self.Ex[i, j] += (dt / (epsilon_0 * self.dy)) * (self.Hz[i, j+1] - self.Hz[i, j])
            self.Ey[i, j] -= (dt / (epsilon_0 * self.dx)) * (self.Hz[i+1, j] - self.Hz[i, j])
            
            
  def _collision(self):
    # TODO Написать обработку коллизий с антенной
    pass
            
            
  def visualize(self, num_steps: int = 100):
    for t in tqdm(range(num_steps)):
      if self.antenna.check_collision:
        self._collision() # недописанная функция
      self._calculate_fields()
      self.Ex[0, 100] += self.source.source_func(t)
      if t % 10 == 0 and t != 0:
        plt.figure(figsize=(6, 6))
        self.antenna.plot_antenna()
        plt.imshow(self.Ex, cmap='hot', extent=[0, self.grid_size * self.dx, 0, self.grid_size * self.dy])
        plt.colorbar(label='Ex (V/m)')
        plt.title(f"Временной шаг {t}")
        plt.xlabel("X (м)")
        plt.ylabel("Y (м)")
        plt.axis("equal")
        plt.tight_layout()
        plt.show()
      
      