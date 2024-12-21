import numpy as np
import matplotlib.pyplot as plt

class Antenna:
  def __init__(self, x0: float, y0: float, diameter: float, focal_length: float):
      self.x0 = x0
      self.y0 = y0
      self.diameter = diameter
      self.focal_length = focal_length
      self.y_sphere = y0 + diameter / 2
      self.x_min = x0 - diameter / 2
      self.x_max = x0 + diameter / 2

    #def check_collision(self, x: float, y: float):
    #    y_parabola = (x - self.x0) ** 2 / (4 * self.focal_length) + self.y0
    #    return np.isclose(y, y_parabola)

  def antenna_equation(self, x: float):
    return self.x_min <= x <= self.x_max and (x - self.x0) ** 2 / (4 * self.focal_length) + self.y0

  def tangent(self, x: float):
    return (x - self.x0) / (2 * self.focal_length)

  def plot_antenna(self, ax):
    xs = np.linspace(self.x_min, self.x_max, 100)
    ys = (xs - self.x0) ** 2 / (4 * self.focal_length) + self.y0

    ax.plot(xs, ys, label="Параболическая антенна", color="blue")
    ax.fill_between(xs, ys, color='lightblue', alpha=0.5)
    ax.scatter(self.x0, self.y0, color="red", label="Центр")
    ax.scatter(self.x0, self.y0 + self.focal_length, color="green", label="Фокус")

    ax.set_title("Параболическая антенна")
    ax.legend()
    ax.grid()