from antenna import Antenna
from source import Source

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from typing import List
import imageio.v2 as imageio

frames = []

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
    def __init__(self, grid_size: int, dx: float, dy: float, antenna: Antenna = None, 
                 sources: List[Source] = [], pml_layers: int = 10, pml_max: float = 1.0):
        self.grid_size = grid_size
        self.dx = dx
        self.dy = dy
        
        # PML layer
        self.total_size = grid_size + 2 * pml_layers
        self.pml_layers = pml_layers
        self.pml_max = pml_max
        
        self.Ex = np.zeros((self.total_size, self.total_size))
        self.Ey = np.zeros((self.total_size, self.total_size))
        self.Hz = np.zeros((self.total_size, self.total_size))
        self.Ex_next = np.zeros((self.total_size, self.total_size))
        self.Ey_next = np.zeros((self.total_size, self.total_size))
        self.Hz_next = np.zeros((self.total_size, self.total_size))
        
        self.sigma_x = np.zeros(self.total_size)
        self.sigma_y = np.zeros(self.total_size)
        
        for i in range(pml_layers):
            x = (pml_layers - i) / pml_layers
            self.sigma_x[i] = self.pml_max * x**2
            self.sigma_x[-(i+1)] = self.pml_max * x**2
                
            self.sigma_y[i] = self.pml_max * x**2
            self.sigma_y[-(i+1)] = self.pml_max * x**2
        
        self.Ex_max = 0
        
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
    
    def renew_focus_field(self):
        self.Ex_max = max(self.Ex_max, self.Ex[int(self.antenna.y0 / self.dy) + int(self.antenna.focal_length / self.dx), int(self.antenna.x0 / self.dx) + self.grid_size//2])
        
    def _calculate_fields(self, dt: float = 0.0025 / (2 * 3e8)):
        for i in range(1, self.total_size):
            for j in range(1, self.total_size):
                damp = 1 / (1 + (self.sigma_x[j] + self.sigma_y[i]) * dt / (2 * EPS_0))
                self.Hz_next[i, j] = damp * (
                    self.Hz_next[i, j] - (dt / MU_0) * (
                        (self.Ey[i, j] - self.Ey[i, j - 1]) / self.dx - 
                        (self.Ex[i, j] - self.Ex[i - 1, j]) / self.dy
                    )
                )
        
        for i in range(self.total_size - 1):
            for j in range(self.total_size - 1):
                x_coord = (j - self.pml_layers - self.total_size//2) * self.dx
                y_coord = (i - self.pml_layers) * self.dy
                
                if self.antenna and self.check_collision(x_coord, y_coord):
                    self.Ex_next[i, j] = -1.0 * self.Ey[i, j] * self.antenna.tangent(x_coord)
                    continue
                    
                if self.antenna and self.check_inside_antenna(x_coord, y_coord):
                    self.Ex_next[i, j] = 0
                    self.Ey_next[i, j] = 0
                    continue
                    
                damp_x = 1 / (1 + self.sigma_x[j] * dt / (2 * EPS_0))
                self.Ex_next[i, j] = damp_x * (
                    self.Ex_next[i, j] + (dt / EPS_0) * (
                        (self.Hz[i + 1, j] - self.Hz[i, j]) / self.dy
                    )
                )
                
                damp_y = 1 / (1 + self.sigma_y[i] * dt / (2 * EPS_0))
                self.Ey_next[i, j] = damp_y * (
                    self.Ey_next[i, j] - (dt / EPS_0) * (
                        (self.Hz[i, j + 1] - self.Hz[i, j]) / self.dx
                    )
                )
        
        self.Ex, self.Ex_next = self.Ex_next, self.Ex
        self.Ey, self.Ey_next = self.Ey_next, self.Ey
        self.Hz, self.Hz_next = self.Hz_next, self.Hz
                
    def visualize(self, num_steps: int = 100, dt: float = 0.0025 / (2 * 3e8), filename: str = "test"):
        source_positions = []
        for source in self.sources: 
            src_i = int(source.y0 / self.dy) + self.pml_layers
            src_j = int(source.x0 / self.dx) + self.total_size//2
            source_positions.append((src_i, src_j))
            self.Ex[src_i, src_j] = source.source_func(0)
            self.Hz[src_i, src_j] = np.sqrt(EPS_0 / MU_0) * source.source_func(0)

        for t in tqdm(range(num_steps)):
            if t < num_steps // 5:
                for (src_i, src_j), source in zip(source_positions, self.sources):
                    self.Ex[src_i, src_j] = source.source_func(t * dt)
                    self.Hz[src_i, src_j] = np.sqrt(EPS_0 / MU_0) * source.source_func(t * dt)
            else:
                self.renew_focus_field()
            
            self._calculate_fields(dt)

            if t % 10 == 0:
                main_Ex = self.Ex[self.pml_layers:-self.pml_layers, self.pml_layers:-self.pml_layers]
                
                fig, ax = plt.subplots(figsize=(6, 6))
                ax.set_xlabel("X (м)")
                ax.set_ylabel("Y (м)")
                ax.set_title("Электромагнитное поле")
                
                if self.antenna:
                    self.antenna.plot_antenna(ax=ax)
                
                min_E = np.min(main_Ex)
                max_E = np.max(main_Ex)
                edge = max(abs(min_E), abs(max_E)) if max(abs(min_E), abs(max_E)) > 0 else 1.0
                
                extent = [
                    -self.grid_size * self.dx / 2,
                    self.grid_size * self.dx / 2,
                    0,
                    self.grid_size * self.dy
                ]
                
                im = ax.imshow(main_Ex, cmap='coolwarm', vmin=-edge, vmax=edge,
                              extent=extent, origin='lower')
                plt.colorbar(im, ax=ax, label='E field strength')
                
                ax.set_title(f"Время {round(t * dt, 10)}")
                ax.axis("equal")
                
                plt.savefig('temp.png')
                plt.close()
                frames.append(imageio.imread('temp.png'))
        
        imageio.mimsave(f'../visuals/{filename}.gif', frames, duration=0.4)