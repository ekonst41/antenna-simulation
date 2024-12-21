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
    
  def _calculate_fields(self, dt: float = 0.0025 / (2 * 3e8)):
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
      fig, ax = plt.subplots(figsize=(6, 6))
      ax.set_xlim(0, self.grid_size * self.dx)  # Установка масштаба по оси X
      ax.set_ylim(0, self.grid_size * self.dy)  # Установка масштаба по оси Y
      ax.set_xlabel("X (м)")
      ax.set_ylabel("Y (м)")
      ax.set_title("Электромагнитное поле")
      ax.axis("equal")

      for t in tqdm(range(num_steps)):
        if self.antenna and self.antenna.check_collision:
          self._collision()  # Обработка коллизий с антенной
        self._calculate_fields()
          
        if (t < num_steps / 10):
          self.Ex[20, 100] += self.source.source_func(t)  # Применение источника

        if t % 10 == 0:
          ax.clear()

          if self.antenna:
            self.antenna.plot_antenna(ax=ax)

                # Отрисовка электрического поля
          ax.imshow(self.Ex, cmap='hot', extent=[self.antenna.x0 - self.grid_size * self.dx / 2, self.antenna.x0 + self.grid_size * self.dx / 2, self.antenna.y0, self.antenna.y0 + self.grid_size * self.dy], origin='lower')
          ax.set_title(f"Временной шаг {t}")
          ax.set_xlabel("X (м)")
          ax.set_ylabel("Y (м)")
          ax.axis("equal")

          # Обновление графика
          plt.draw()
          plt.pause(0.01)  # Пауза для обновления

      plt.show()
      
      